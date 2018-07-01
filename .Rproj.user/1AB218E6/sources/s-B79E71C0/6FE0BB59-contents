################################################################################
# Complex-valued EMVS: Simulation study 2                                      #
# The code is based on general complex normal and embed C++ using Rcpp         #
# Cheng-Han Yu, UC Santa Cruz                                                  #
################################################################################
rm(list = ls())
################################################################################
# load required packages
################################################################################
library(fields)
library(Matrix) # for sparse matrix (sparseMatrix)
library(doMC) # simple Parallel computing
library(mvnfast)
library(abind)
library(Hmisc)

################################################################################
# source required files
################################################################################
source("EMVSmainfcn.R")
source("./CVEMVS/log_g_sig_noncir.R")
source("./CVEMVS/log_prior.R") # same as circular
source("./CVEMVS/E_beta_binom_noncir.R")   # E-step
source("./CVEMVS/M_beta_noncir.R") # M-step
source("./CVEMVS/M_sigma_sig_noncir.R") # M-step
source("./CVEMVS/M_theta.R")  # M-step  same as circular

source("./CVEMVS/EMVSmainfcn_noncir1.R")
source("./CVEMVS/EMVSbest.R") # same as circular
source("./CVEMVS/EMVSplot.R") # same as circular

source("./CVEMVS/CEMVSfcnscpp.R")
source("./data/loadData.R")

################################################################################
# Run the general CV EMVS algorithms 
################################################################################
# set parameters
################
q <- 1
nn <- n.time * 2
qq <- q * 2
# v0s <- 0.006
v0s <- 0.001 * 1:15
v1 <- 1
epsilon <- 1e-3
thre <- c(0.5, 0.8722)

# CV algorithm
################
system.time(sim2.em.cpp <- EM_spikecir_cpp_full(X.cplx, Y.cv, v0s, v1 = 1, epsilon = epsilon,
                                                a_th = 1, b_th = 1, thre = thre[1], 
                                                spike.cor = 0, slab.cor = 0))

# be careful about the index (q, Nvoxels, L)
# pprob <- sim2.em.cpp$pprob[1, , which(v0s == EMVSbest(sim2.em.cpp)$v0)]




# Plotting
################
par(mfrow = c(1, 3), omi = c(0, 0, 0, 0), mar = c(.1, .1, 1.5, .1))
data.mean.list <- lapply(list('Re' = dataR.r, 'Im' = dataI.r, 'Mod' = dataM.r, 
                              'Arg' = dataA.r), function(x) apply(x, c(2, 3), mean))
# Activation maps
img_fmri(output = sim2.em.cpp, thre = 0.5, v0s = v0s, v1 = v1, iscplx = TRUE, 
         isstr = FALSE, cex.pt = 1, cpp = TRUE)
# Strength map
img_fmri(output = sim2.em.cpp, thre = 0.5, v0s = v0s, v1 = v1, iscplx = TRUE, 
         isstr = TRUE, title = "CV Estimated Strength", isbar = FALSE,
         cex.main = 1.5, cex.pt = 1, cpp = TRUE)

# Strength map with the color bar
img_fmri(output = sim2.em.cpp, thre = 0.5, v0s = v0s, v1 = v1, iscplx = TRUE, 
         isstr = TRUE, title = "CV Estimated Strength", isbar = TRUE,
         cex.main = 1.5, cex.pt = 1, cpp = TRUE)














