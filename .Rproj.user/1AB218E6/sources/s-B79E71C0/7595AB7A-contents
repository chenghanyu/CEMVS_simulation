################################################################################
# load and shape the data
################################################################################
# MATLAB Code to save the data in txt format that can be read into R
# ==========
# data = load('GrantPrelimSimData.mat')
# datafield = fieldnames(data);
# dlmwrite('simulation2.txt', data.(datafield{1}));
# ==========

# read the data into R and save it as simulation2.RData
# simdata <- read.table("./data/simulation2.txt", sep = ",")
# save(simdata, file = "./data/simulation2.RData")

load("./data/simulation2.RData")  # load the simulated data object called simdata
N <- 96
n.time <- 490
data <- as.matrix(simdata)
dataR <- array(Re(data), dim = c(N, N, n.time))
dataI <- array(Im(data), dim = c(N, N, n.time))
dataR <- aperm(dataR, c(3, 1, 2))  # array permutation # 490 * 96 * 96
dataI <- aperm(dataI, c(3, 1, 2))
dataC <- dataR + 1i * dataI
dataM <- Mod(dataC)
dataA <- Arg(dataC)

dataM.r <- dataA.r <- dataR.r <- dataI.r <- array(0, c(n.time, N, N))
for (t in 1:n.time) {
    dataM.r[t, , ] <- rotate(dataM[t, , ])
    dataA.r[t, , ] <- rotate(dataA[t, , ])
    dataR.r[t, , ] <- rotate(dataR[t, , ])
    dataI.r[t, , ] <- rotate(dataI[t, , ])
}


################################################################################
# generate the covariate and reference function with time trend 
################################################################################
library(neuRosim)
n <- 510
# get rid of the first 20 time points
ones <- rep(1, n.time) # intercept
x1 <- 1:n.time  # trend
totaltime <- n
onsets <- seq(30, n, by = 30)
dur <- 15
s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, 
                  accuracy = 1)
# HRF
canonical <- specifydesign(totaltime = totaltime, onsets = list(onsets),
                           durations = list(dur), effectsize = 1, TR = 1, 
                           conv = "double-gamma")
# design matrix
x <- cbind(ones, x1 / n.time, canonical[21:n])
X <- as.matrix(bdiag(x, x))
x.re <- canonical[21:n]
# first dimension is time
# combine real and imaginary part by time
Y.cv <- abind(dataR.r, dataI.r, along = 1)  

# remove intercept and trend
xstar <- cbind(ones, x1)
Xstar <- X[, c(1, 2, 4, 5)]
IIc <- (diag(2 * n.time) - Xstar %*% solve(t(Xstar) %*% Xstar) %*% t(Xstar))
II <- (diag(n.time) - xstar %*% solve(t(xstar) %*% xstar) %*% t(xstar))
X.cplx <- (diag(2 * n.time) - Xstar %*% solve(t(Xstar) %*% Xstar) %*% t(Xstar)) %*% 
    X[, c(3, 6)]
x.re <- (diag(n.time) - xstar %*% solve(t(xstar) %*% xstar) %*% t(xstar)) %*% x.re
x.mod <- x.re
Y.mod <- dataM.r

Y.cv <- apply(Y.cv, c(2, 3), function(x) IIc %*% x)
Y.mod <- apply(Y.mod, c(2, 3), function(x) II %*% x)
