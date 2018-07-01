################################################################################
# Complex-valued fMRI Simulation Study: Non-circular Analysis                  #
# EMVS one theta and common variance                                           #
# No trend data                                                                #
# Cheng-Han Yu, UC Santa Cruz                                                  #
################################################################################
rm(list = ls())
# load required packages
################################################################################
library(fields)
library(Matrix) # for sparse matrix (sparseMatrix)
library(doMC) # simple Parallel computing
library(mvnfast)
library(abind)
# source required files
################################################################################
source("log_g_sig_noncir.R")
source("E_beta_binom_noncir.R")
source("log_prior.R") # same as circular
source("M_beta_noncir.R")
source("M_theta.R")  # same as circular
source("M_sigma_sig_noncir.R")
source("EMVSbest.R") # same as circular
source("EMVSplot.R") # same as circular
source("tpfpfcn.R") # same as circular
source("EMVSmainfcn_noncir.R")
source("EMVSmainfcn.R")
source("log_g_sig.R")
source("E_beta_binom.R")
source("M_beta.R")
source("M_sigma_sig.R")

# Experimental Design (Same as circular model)
################################################################################
source("./Data/SimExpDesign.R")

# Now the script below only containa functions for Generating Design matrix and 
# noncircular coefficients
################################################################################
source("./Data/SimDesignMatNoncirCoef.R") # now only contain X and functions

# Make it simple. Consider one correlation and one set of regression coef mean 
# at each time only.
# 1. ignore sig2 when generating reg. coeff.

# Generate Regression coefficients
########################
phi <- pi / 4 
(alpha1 <- cos(phi))
(alpha2 <- sin(phi))
sig2 <- 0.25
SNR <- c(0.5)
beta0 <- sqrt(sig2) * SNR
cor.ri <- 0  # correlation for active voxels
cor.ri.nonact <- 0  # correlation for nonactive voxels

# Create covriance matrix of reg coef Gamma using function create_covmat_ga()
# =========================================
Vrr <- 1  # When Vrr = Vii = 1, covriance matrix is equivalent to correlation matrix
Vii <- 1
Vrr_nonact <- 1
Vii_nonact <- 1
covmat <- create_covmat_ga(Vrr, Vii, Vrr_nonact, Vii_nonact, 
                           cor.ri = cor.ri, cor.ri.nonact, 
                           alpha1, alpha2, p = 1, sig2 = sig2, issig2 = TRUE)
(Sig_ga <- covmat$Sig_ga)
(Sig_ga_nonact <- covmat$Sig_ga_nonact)

set.seed(7826)
beta.mu.re <- 1  # length one for simplicity
beta.mu.im <- 1  # could be different
(ga.mu.re <- alpha1 * beta.mu.re)
(ga.mu.im <- alpha2 * beta.mu.im)

# Generate reg. coeff using function create_regcoef()
# =========================================

# This generates coefficients for active and nonactive voxels. 
# One at each voxel and coefficient values are correlated based on Sig_ga
gacoef <- create_ga_bold(ga.mu.re = ga.mu.re, 
                         ga.mu.im = ga.mu.im, 
                         Sig_ga, Sig_ga_nonact,
                         final, alpha1, alpha2)
par(mfrow = c(2,1))
plot(gacoef$ga[, 1], gacoef$ga[, 2], xlab = "ga.re", ylab = "ga.im",
     main = paste("Act ga.mu.re =", round(ga.mu.re[1],2),  
                  "ga.mu.im =", round(ga.mu.im[1],2),  
                  "; cor =", cor.ri))

plot(gacoef$ga_nonact[, 1], gacoef$ga_nonact[, 2], 
     xlab = "ga.re", ylab = "ga.im",
     main = paste("NonAct ga.mu.re = 0 ga.mu.im = 0",  
                  "; cor =", cor.ri.nonact))


# Generate reg. coeff vector including the baseline 
# using function create_regcoef_vec()
# =========================================

# Here we actually force coefficients of nonactive voxels to be exactly zero
ga <- gacoef$ga
ga_nonact <- gacoef$ga_nonact
Ga.vec <- create_regcoef_vec(ga, ga_nonact, final, beta0, alpha1, alpha2)

Ga.re <- Ga.vec$Ga.re
Ga.im <- Ga.vec$Ga.im


# Simulate data sets
################################################################################
source("./Data/GenNonCirData1.R") # p = 1



# set tunning parameters
########################
q <- ncol(x.mod)
v0s <- v0s.mod <- seq(0.001, 0.03, length = 10)
v1 <- 1
thre <- c(0.5, 0.8722, 0.95)
epsilon <- 1e-3



# Run the algorithm: function for different cor.ri
# ==============================
cor.ri.mdl <- 0.9

em_theta_sig_noncir_cor <- EM_theta_sig_noncir(X.cplx, Y, v0s, v1 = 1, 
                                               epsilon = epsilon,
                                               a = 1, b = 1, 
                                               thre = thre[1], 
                                               cor.ri = cor.ri.mdl)

em_theta_sig_cir_cor <- EM_theta_sig(X.cplx, Y, v0s, v1 = 1, epsilon = epsilon,
                                     a = 1, b = 1, thre = thre[1])


em_theta_sig_spikecir_cor <- EM_theta_sig_spikecir(X.cplx, Y, v0s, v1 = 1, 
                                                   epsilon = epsilon,
                                                   a = 1, b = 1, thre = thre[1], 
                                                   cor.ri = cor.ri.mdl)

save(em_theta_sig_noncir_cor, em_theta_sig_cir_cor, em_theta_sig_spikecir_cor, 
     file = "em_res_cor.RData")
# load("em_res_cor.RData", verbose = TRUE)

################################################################################
imageplot <- function(aa, posprob, main) {
    for (i in aa) {
        image(posprob > i, main = paste(main), axes = FALSE,
              col = c("navy", "red3"), cex.main = 1.5)
    }
}

coef_plot <- function(Gamma.re, Gamma.im, res_noncir = NULL, res_cir = NULL, 
                      res_spikecir = NULL, final,
                      embest_non = NULL, embest = NULL, embest_spike = NULL, 
                      cor.ri, cex = 0.45, main = "") {
    # par(mfrow = c(2, 5))
    par(mar = c(4,4,2,1))
    xlim <- c(range(Gamma.re)[1] - 0.3, range(Gamma.re)[2] + 0.2)
    ylim <- c(range(Gamma.im)[1] - 0.3, range(Gamma.im)[2] + 0.2)
    plot(Gamma.re, Gamma.im, xlim = xlim, ylim =ylim, 
         main = paste("true coef; cor.ri =", cor.ri, main), pch = 19, cex = cex, col = ifelse(final != 0, "red", "black"))
    legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
    if(!is.null(res_noncir)) {
        plot(res_noncir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_noncir$betas[1, 2, , ], col = ifelse(final != 0, "red", "black"),
             cex = cex, xlim = xlim, ylim = ylim, main = "noncir coeff", pch = 19)
        legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
        Beta_init <- res_noncir$Beta_init
    }
    
    if(!is.null(res_cir)) {
        plot(res_cir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_cir$betas[1, 2, , ], col = ifelse(final != 0, "red", "black"),
             cex = cex, xlim = xlim, ylim = ylim, main = "cir coeff", pch = 19)
        legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
        Beta_init <- res_cir$Beta_init
    }
    
    if (!is.null(res_spikecir)) {
        plot(res_spikecir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_spikecir$betas[1, 2, , ], col = ifelse(final != 0, "red", "black"),
             cex = cex, xlim =xlim, ylim = ylim, main = "spikecir coeff", pch = 19)
        legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
        Beta_init <- res_spikecir$Beta_init
    }
    
    
    plot(Beta_init[1, , ], Beta_init[2, , ], cex = cex,
         xlim = xlim, ylim = ylim, main = "lse", pch = 19,
         xlab = "re", ylab = "im", col = ifelse(final != 0, "red", "black"))
    legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
    
    
    plot(Gamma.re, Gamma.im, xlim =xlim, ylim = ylim, 
         main = "true coef", pch = 19, cex = cex, col = ifelse(final != 0, "red", "black"))
    legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
    
    if(!is.null(res_noncir)) {
        colmat<- matrix("black", N, N)
        # blue: true positives
        # purple: false positives
        # green: false negatives
        colmat[res_noncir$pprob[embest_non, 1, , ] > 0.5 & final != 0] <- "blue"
        colmat[res_noncir$pprob[embest_non, 1, , ] > 0.5 & final == 0] <- "purple"
        colmat[res_noncir$pprob[embest_non, 1, , ] < 0.5 & final != 0] <- "green"
        plot(res_noncir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_noncir$betas[1, 2, , ], 
             col = colmat,
             cex = cex, xlim = xlim, ylim = ylim, main = "noncir coeff label", pch = 19)
        legend("bottomright", c("true-postive", "false-positive", "false negative"), 
               col = c("blue", "purple", "green"), bty = "n", pch = 19)
    }
    
    if(!is.null(res_cir)) {
        colmat<- matrix("black", N, N)
        colmat[res_cir$pprob[embest, 1, , ] > 0.5 & final != 0] <- "blue"
        colmat[res_cir$pprob[embest, 1, , ] > 0.5 & final == 0] <- "purple"
        colmat[res_cir$pprob[embest, 1, , ] < 0.5 & final != 0] <- "green"
        plot(res_cir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_cir$betas[1, 2, , ], 
             col = colmat,
             cex = cex, xlim = xlim, ylim = ylim, main = "cir coeff label", pch = 19)
        legend("bottomright", c("true-postive", "false-positive", "false negative"), 
               col = c("blue", "purple", "green"), bty = "n", pch = 19)
    }
    
    if (!is.null(res_spikecir)) {
        colmat<- matrix("black", N, N)
        colmat[res_spikecir$pprob[embest_spike, 1, , ] > 0.5 & final != 0] <- "blue"
        colmat[res_spikecir$pprob[embest_spike, 1, , ] > 0.5 & final == 0] <- "purple"
        colmat[res_spikecir$pprob[embest_spike, 1, , ] < 0.5 & final != 0] <- "green"
        plot(res_spikecir$betas[1, 1, , ], xlab = "re", ylab = "im",
             res_spikecir$betas[1, 2, , ], 
             col = colmat,
             cex = cex, xlim =xlim, ylim = ylim, main = "spikecir coeff label", pch = 19)
        legend("bottomright", c("true-postive", "false-positive", "false negative"), 
               col = c("blue", "purple", "green"), bty = "n", pch = 19)
    }
    
    plot(Beta_init[1, , ], Beta_init[2, , ], cex = cex,
         xlim = xlim, ylim = ylim, main = "lse", pch = 19,
         xlab = "re", ylab = "im", col = ifelse(final != 0, "red", "black"))
    legend("bottomright", c("active", "non-active"), col = c("red", "black"), bty = "n", pch = 19)
}


################################################################################

# ===================
# activation plot
# --------------
par(mfrow = c(1, 3))
par(mar = c(1,1,2,1))
embest_non.vec <- which(v0s == EMVSbest(em_theta_sig_noncir_cor)$v0)
imageplot(thre[1], em_theta_sig_noncir_cor$pprob[embest_non.vec, 1, , ], 
          main = paste("CV-noncir cor=", cor.ri))
embest_cir.vec <- which(v0s == EMVSbest(em_theta_sig_cir_cor)$v0)
imageplot(thre[1], em_theta_sig_cir_cor$pprob[embest_cir.vec, 1, , ], 
          main = paste("CV-cir cor=", cor.ri))
embest_spike.vec <- which(v0s == EMVSbest(em_theta_sig_spikecir_cor)$v0)
imageplot(thre[1], em_theta_sig_spikecir_cor$pprob[embest_spike.vec, 1, , ], 
          main = paste("CV-spikecir cor=", cor.ri))

# coefficient plot
# --------------
par(mfrow = c(2, 5))
coef_plot(Gamma.re = Ga.re[-1], 
          Gamma.im = Ga.im[-1], 
          res_noncir = em_theta_sig_noncir_cor, 
          res_cir = em_theta_sig_cir_cor, 
          res_spikecir = em_theta_sig_spikecir_cor, 
          final = final,
          embest_non = embest_non.vec, 
          embest = embest_cir.vec, 
          embest_spike = embest_spike.vec, 
          cor.ri = cor.ri)



# MSE of beta
# ------------
compute_mse <- function(res.non, res.cir, res.spike, 
                        embest.non, embest.cir, embest.spike, 
                        cor.ri, Ga.re, Ga.im) {
    MSE.cor.non <- matrix(0, 1, 3)
    rownames(MSE.cor.non) <- paste0("cor", cor.ri)
    colnames(MSE.cor.non) <- c("all", "true_act", "label_act")
    MSE.cor.cir <- matrix(0, 1, 3)
    rownames(MSE.cor.cir) <- paste0("cor", cor.ri)
    colnames(MSE.cor.cir) <- c("all", "true_act", "label_act")
    MSE.cor.spike <- matrix(0, 1, 3)
    rownames(MSE.cor.spike) <- paste0("cor", cor.ri)
    colnames(MSE.cor.spike) <- c("all", "true_act", "label_act")
    # =============
    betas <- res.non$betas[embest.non, , , ]
    # all
    MSE.cor.non[1, 1] <- mean((as.vector(betas[1, , ]) - Ga.re[-1]) ^ 2 + 
                                  (as.vector(betas[2, , ]) - Ga.im[-1]) ^ 2) 
    # true act
    MSE.cor.non[1, 2] <- mean((as.vector(betas[1, , ])[final != 0] - Ga.re[-1][final != 0]) ^ 2 + 
                                  (as.vector(betas[2, , ])[final != 0] - Ga.im[-1][final != 0]) ^ 2) 
    # labeled act
    labeled <- res.non$pprob[embest.non, 1, , ] > 0.5
    MSE.cor.non[1, 3] <- mean((as.vector(betas[1, , ])[labeled] - Ga.re[-1][labeled]) ^ 2 + 
                                  (as.vector(betas[2, , ])[labeled] - Ga.im[-1][labeled]) ^ 2) 
    # =============
    betas <- res.cir$betas[embest.cir, , , ]
    # all
    MSE.cor.cir[1, 1] <- mean((as.vector(betas[1, , ]) - Ga.re[-1]) ^ 2 + 
                                  (as.vector(betas[2, , ]) - Ga.im[-1]) ^ 2) 
    # true act
    MSE.cor.cir[1, 2] <- mean((as.vector(betas[1, , ])[final != 0] - Ga.re[-1][final != 0]) ^ 2 + 
                                  (as.vector(betas[2, , ])[final != 0] - Ga.im[-1][final != 0]) ^ 2) 
    # labeled act
    labeled <- res.cir$pprob[embest.cir, 1, , ] > 0.5
    MSE.cor.cir[1, 3] <- mean((as.vector(betas[1, , ])[labeled] - Ga.re[-1][labeled]) ^ 2 + 
                                  (as.vector(betas[2, , ])[labeled] - Ga.im[-1][labeled]) ^ 2) 
    # =============
    betas <- res.spike$betas[embest.spike, , , ]
    # all
    MSE.cor.spike[1, 1] <- mean((as.vector(betas[1, , ]) - Ga.re[-1]) ^ 2 + 
                                  (as.vector(betas[2, , ]) - Ga.im[-1]) ^ 2) 
    # true act
    MSE.cor.spike[1, 2] <- mean((as.vector(betas[1, , ])[final != 0] - Ga.re[-1][final != 0]) ^ 2 + 
                                  (as.vector(betas[2, , ])[final != 0] - Ga.im[-1][final != 0]) ^ 2) 
    # labeled act
    labeled <- res.spike$pprob[embest.spike, 1, , ] > 0.5
    MSE.cor.spike[1, 3] <- mean((as.vector(betas[1, , ])[labeled] - Ga.re[-1][labeled]) ^ 2 + 
                                  (as.vector(betas[2, , ])[labeled] - Ga.im[-1][labeled]) ^ 2) 
    
    return(list(MSE.cor.non = MSE.cor.non,
                MSE.cor.cir = MSE.cor.cir,
                MSE.cor.spike = MSE.cor.spike))
}

MSEres <- compute_mse(res.non = em_theta_sig_noncir_cor, 
            res.cir = em_theta_sig_cir_cor, 
            res.spike = em_theta_sig_spikecir_cor, 
            embest.non = embest_non.vec, 
            embest.cir = embest_cir.vec, 
            embest.spike = embest_spike.vec, 
            cor.ri = cor.ri, 
            Ga.re = Ga.re, 
            Ga.im = Ga.im)

MSEres$MSE.cor.non
MSEres$MSE.cor.cir
MSEres$MSE.cor.spike

# ################################################################################
# # different beta
# # ===================
# # activation plot
# # ------------
# par(mfrow = c(3, 3))
# for (s in 1:3) {
#     embest_non.vec[s] <- which(v0s == EMVSbest(em_theta_sig_noncir_coef[[s]])$v0)
#     imageplot(thre[1], em_theta_sig_noncir_coef[[s]]$pprob[embest_non.vec[s], 1, , ], 
#               main = paste("CV-noncir beta=", beta.mu.vec[s]))
#     # image.plot(em_theta_sig_noncir_cor[[s]]$pprob[embest, 1, , ], 
#     #            main = "CV-noncir", axes = FALSE, cex.main = 2)
#     embest_cir.vec[s] <- which(v0s == EMVSbest(em_theta_sig_cir_coef[[s]])$v0)
#     imageplot(thre[1], em_theta_sig_cir_coef[[s]]$pprob[embest_cir.vec[s], 1, , ], 
#               main = paste("CV-cir beta=", beta.mu.vec[s]))
#     embest_spike.vec[s] <- which(v0s == EMVSbest(em_theta_sig_spikecir_coef[[s]])$v0)
#     imageplot(thre[1], em_theta_sig_spikecir_coef[[s]]$pprob[embest_spike.vec[s], 1, , ], 
#               main = paste("CV-spikecir beta=", beta.mu.vec[s]))
# }
# 
# # coefficient plot
# # ------------
# for (s in 1:3) {
#     coef_plot(Gamma.re = coef.mu.cor.lst[[s]][[3]]$Beta.re[-1], 
#               Gamma.im = coef.mu.cor.lst[[s]][[3]]$Beta.im[-1], 
#               res_noncir = em_theta_sig_noncir_coef[[s]], 
#               res_cir = em_theta_sig_cir_coef[[s]], 
#               res_spikecir = em_theta_sig_spikecir_coef[[s]], 
#               final = final,
#               embest_non = embest_non.vec[s], 
#               embest = embest_cir.vec[s], 
#               embest_spike = embest_spike.vec[s], 
#               cor.ri = cor.ri.vec[3],
#               cex = 0.25)
# }
# 
# # MSE of beta
# # ------------
# MSE.coef.non <- matrix(0, 3, 3)
# rownames(MSE.coef.non) <- paste0("coef", beta.mu.vec)
# colnames(MSE.coef.non) <- c("all", "true_act", "label_act")
# 
# MSE.coef.cir <- matrix(0, 3, 3)
# rownames(MSE.coef.cir) <- paste0("cor", cor.ri.vec)
# colnames(MSE.coef.cir) <- c("all", "true_act", "label_act")
# 
# MSE.coef.spike <- matrix(0, 3, 3)
# rownames(MSE.coef.spike) <- paste0("cor", cor.ri.vec)
# colnames(MSE.coef.spike) <- c("all", "true_act", "label_act")
# 
# for (s in 1:3) {
#     betas <- em_theta_sig_noncir_coef[[s]]$betas[embest_non.vec[s], , , ]
#     # all
#     MSE.coef.non[s, 1] <- mean((as.vector(betas[1, , ]) - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1]) ^ 2 + 
#                                    (as.vector(betas[2, , ]) - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1]) ^ 2) 
#     # true act
#     MSE.coef.non[s, 2] <- mean((as.vector(betas[1, , ])[final != 0] - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][final != 0]) ^ 2 + 
#                                    (as.vector(betas[2, , ])[final != 0] - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][final != 0]) ^ 2) 
#     # labeled act
#     labeled <- em_theta_sig_noncir_coef[[s]]$pprob[embest_non.vec[s], 1, , ] > 0.5
#     MSE.coef.non[s, 3] <- mean((as.vector(betas[1, , ])[labeled] - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][labeled]) ^ 2 + 
#                                    (as.vector(betas[2, , ])[labeled] - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][labeled]) ^ 2) 
#     ##################
#     betas <- em_theta_sig_cir_coef[[s]]$betas[embest_cir.vec[s], , , ]
#     # all
#     MSE.coef.cir[s, 1] <- mean((as.vector(betas[1, , ]) - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1]) ^ 2 + 
#                                    (as.vector(betas[2, , ]) - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1]) ^ 2) 
#     # true act
#     MSE.coef.cir[s, 2] <- mean((as.vector(betas[1, , ])[final != 0] - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][final != 0]) ^ 2 + 
#                                    (as.vector(betas[2, , ])[final != 0] - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][final != 0]) ^ 2) 
#     # labeled act
#     labeled <- em_theta_sig_cir_coef[[s]]$pprob[embest_cir.vec[s], 1, , ] > 0.5
#     MSE.coef.cir[s, 3] <- mean((as.vector(betas[1, , ])[labeled] - 
#                                     coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][labeled]) ^ 2 + 
#                                    (as.vector(betas[2, , ])[labeled] - 
#                                         coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][labeled]) ^ 2)
#     ##################
#     betas <- em_theta_sig_spikecir_coef[[s]]$betas[embest_spike.vec[s], , , ]
#     # all
#     MSE.coef.spike[s, 1] <- mean((as.vector(betas[1, , ]) - 
#                                       coef.mu.cor.lst[[s]][[3]]$Beta.re[-1]) ^ 2 + 
#                                      (as.vector(betas[2, , ]) - 
#                                           coef.mu.cor.lst[[s]][[3]]$Beta.im[-1]) ^ 2) 
#     # true act
#     MSE.coef.spike[s, 2] <- mean((as.vector(betas[1, , ])[final != 0] - 
#                                       coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][final != 0]) ^ 2 + 
#                                      (as.vector(betas[2, , ])[final != 0] - 
#                                           coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][final != 0]) ^ 2) 
#     # labeled act
#     labeled <- em_theta_sig_spikecir_coef[[s]]$pprob[embest_spike.vec[s], 1, , ] > 0.5
#     MSE.coef.spike[s, 3] <- mean((as.vector(betas[1, , ])[labeled] - 
#                                       coef.mu.cor.lst[[s]][[3]]$Beta.re[-1][labeled]) ^ 2 + 
#                                      (as.vector(betas[2, , ])[labeled] - 
#                                           coef.mu.cor.lst[[s]][[3]]$Beta.im[-1][labeled]) ^ 2)
# }







