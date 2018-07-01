################################################################################
# Complex-valued EMVS: Simulation study 2                                      #
# The code is based on circular normal and is in PURE R environment            #                         
# Cheng-Han Yu, UC Santa Cruz                                                  #
################################################################################
rm(list = ls())
################################################################################
# Load associated required packages
################################################################################
library(Matrix)
library(mvtnorm) 
library(neuRosim)
library(abind)
library(fields)
library(Hmisc)
library(memoise)
library(R.matlab)

################################################################################
# Source files for the algorithm functions and loading data (Y, X)
################################################################################
# functions
source("EMVSmainfcn.R")
source("./circularEMVS/M_beta.R")
source("./circularEMVS/log_prior.R")
source("./circularEMVS/log_g_sig.R")
source("./circularEMVS/E_beta_binom.R")
source("./circularEMVS/M_theta.R")
source("./circularEMVS/M_sigma_sig.R")
source("./circularEMVS/EMVSbest.R")
# load data (Y, X) 
source("./data/loadData.R")

################################################################################
# Run MO and CV EM algorithms 
################################################################################
# set parameters
################
q <- 1
nn <- n.time * 2
qq <- q * 2
v0s <- 0.006  # the optimal v0 chosen from a grid
# v0s <- 0.001 * 1:15
v1 <- 1
epsilon <- 1e-3

# CV
###############
# EM_theta_sig_test <- memoise(EM_theta_sig)

sim2.em <- EM_theta_sig(X.cplx, Y.cv, v0s, v1, epsilon = epsilon, 
                        a = 1, b = 1, thre = 0.5)
# MO
##########
sim2.em.mod <- EMmod_theta_sig(x.mod, Y.mod, v0s, v1, epsilon = epsilon, 
                               a = 1, b = 1, thre = 0.5)


################################################################################
# Reproduce the same Figure 9 and Figure 10 in the manuscript
################################################################################
#-------------------------------------------------------------------------------
# Figure 9: 
# Left: True activation map. 
# Middle: Activation map from complex-valued EM at optimal $v_0$. 
# Right: Activation map from magnitude-only EM at optimal $v_0$
#-------------------------------------------------------------------------------
png("./figures/sim2_act.png", width = 1000, height = 350)
par(mfrow = c(1, 3), omi = c(0, 0, 0, 0), mar = c(.1, .1, 1.5, .1))
data.mean.list <- lapply(list('Re' = dataR.r, 'Im' = dataI.r, 'Mod' = dataM.r, 
                             'Arg' = dataA.r), function(x) apply(x, c(2, 3), mean))
image(data.mean.list$Mod, col = gray((0:N) / N), axes = FALSE,
      main = "True activation map", cex.main = 1.2)

# two activation regions
cex.pt <- 1.1
points(matrix(rep(43:47, 5) - 0.5, 5, 5) / N, 
       matrix(rep(43:47, each = 5)- 0.5, 5, 5) / N, col = "red", pch = 15, 
       cex = cex.pt)
points(matrix(rep(31:35, 5)- 0.5, 5, 5) / N, 
       matrix(rep(38:42, each = 5)- 0.5, 5, 5) / N, col = "red", pch = 15, 
       cex = cex.pt)

# Activation maps
img_fmri(output = sim2.em, thre = 0.5, v0s = v0s, v1 = v1, iscplx = TRUE, 
         isstr = FALSE, cex.pt = cex.pt)
img_fmri(output = sim2.em.mod, thre = 0.5, v0s = v0s, v1 = v1, iscplx = FALSE, 
         isstr = FALSE, cex.pt = cex.pt)

dev.off()

#-------------------------------------------------------------------------------
# Figure 10: (with margin)
# Left: True strength map. 
# Middle: Strength map from complex-valued EM at optimal $v_0$. 
# Right: Strength map from magnitude-only EM at optimal $v_0$
#-------------------------------------------------------------------------------
png("./figures/sim2brainstr.png", width = 1000, height = 350)
par(mfrow = c(1, 3), omi = c(0, 0, 0, 0), mar = c(.1, .1, 1.5, .1))
layout(matrix(c(rep(1, 5), rep(2, 5), rep(3, 6)), nrow = 1, ncol = 16, 
              byrow = TRUE))

# magnitude contrast strength
beta1s <- readMat('./data/beta1s.mat')
beta1s.mat <- rotate(matrix(beta1s$beta1s, N, N))
true.strength.scale5 <- scale_fcn(beta1s.mat, 0, 0.5)

image(data.mean.list$Mod, col = gray((0:N) / N), axes = FALSE, 
      main = "True Strength", cex.main = 2)
addTrueVoxelStr(true.strength.scale5)
cex.pt <- 1
img_fmri(output = sim2.em, thre = 0.5, v0s = v0s, v1 = v1, iscplx = TRUE, 
         isstr = TRUE, title = "CV Estimated Strength", isbar = FALSE,
         cex.main = 2, cex.pt = cex.pt)
# par(mar = c(.85, .1, 1.2, 5))
img_fmri(output = sim2.em.mod, thre = 0.5, v0s = v0s, v1 = v1, iscplx = FALSE, 
         isstr = TRUE, title = "MO Estimated Strength", isbar = TRUE,
         bar.size = c(0.1, 4.5), cex.main = 2, cex.pt = cex.pt)
dev.off()

#-------------------------------------------------------------------------------
# Figure 10: (without margin)
# Left: True activation map. 
# Middle: Activation map from complex-valued EM at optimal $v_0$. 
# Right: Activation map from magnitude-only EM at optimal $v_0$
#-------------------------------------------------------------------------------
png("./figures/sim2brainstr_cut.png", width = 1000, height = 350)
# Cut margins
par(mfrow = c(1, 3), omi = c(0, 0, 0, 0), mar = c(.1, .1, 1.5, .1))
layout(matrix(c(rep(1, 5), rep(2, 5), rep(3, 6)), nrow = 1, ncol = 16, byrow = TRUE))
cut.row.idx <- c(1:20, 77:96)
cut.col.idx <- c(1:15, 82:96)
img.cut <- data.mean.list$Mod[-cut.row.idx, -cut.col.idx] # dim = c(56, 66)
image(img.cut, col = gray((0:96) / 96), axes = FALSE, 
      main = "True Strength", cex.main = 2)
true.strength.scale5.cut <- true.strength.scale5[-cut.row.idx, -cut.col.idx]
cex.pt <- 1.5
addTrueVoxelStr(true.strength.scale5.cut, cex.pt = cex.pt)

img_fmri(braindata = img.cut, output = sim2.em, thre = 0.5, v0s = v0s, v1 = v1, 
         iscplx = TRUE, isstr = TRUE, title = "CV Estimated Strength", 
         isbar = FALSE, iscut = TRUE, bar.size = c(0.1, 4.5),
         cex.pt = cex.pt, cex.main = 2)
img_fmri(braindata = img.cut, output = sim2.em.mod, thre = 0.5, v0s = v0s, v1 = v1, 
         iscplx = FALSE, isstr = TRUE, title = "MO Estimated Strength",
         iscut = TRUE, cex.pt = cex.pt, bar.size = c(0.1, 4.5),
         isbar = TRUE, cex.main = 2)
dev.off()




