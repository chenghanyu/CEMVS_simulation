###############################################################################
# general CEMVS functions using Rcpp
###############################################################################
library(compiler)
library(Rcpp)
library(RcppArmadillo)
library(inline)

# souce cpp files
sourceCpp("./CVEMVS/rcpp/E_beta_binom_spikecir_cpp.cpp")
sourceCpp("./CVEMVS/rcpp/M_beta_noncir_cpp.cpp")
sourceCpp("./CVEMVS/rcpp/M_theta_cpp.cpp")
sourceCpp("./CVEMVS/rcpp/M_sigma_sig_noncir_cpp.cpp")

# for full cpp
sourceCpp("./CVEMVS/rcpp/EM_spikecir_main_cpp.cpp")
sourceCpp("./CVEMVS/rcpp/log_g_sig_noncir_cpp.cpp")
sourceCpp("./CVEMVS/rcpp/log_prior_cpp.cpp")

EM_spikecir_cpp <- function(X, Y, v0s, v1, epsilon, theta, a_th, b_th, thre, 
                            spike.cor, slab.cor) {
    # X: n x qq matrix 
    # Y: n x N x N Array
    # v0: vector of length L
    # v1: scalar
    # beta_init: qq x 1 vector
    # sigma_init: scalar
    # epsilon: scalar
    # theta: scalar
    # a: scalar
    # b: scalar
    # va_g: scalar
    #################################################################
    if (missing(theta)){theta = 0.5}
    
    if (class(X) != "matrix") {
        tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("X must be a matrix or able to be coerced to a matrix")
    }
    
    if (class(Y) != "array") {
        tmp <- try(Y <- as.array(Y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("Y must be an array or able to be coerced to an array")
    }
    
    if (sum(v1<v0s) > 0)
        stop("v1 has to be larger than v0")
    
    if (length(a_th) > 1) 
        stop("a has to be a scalar")
    
    if (length(b_th) > 1) 
        stop("b has to be a scalar")
    
    if (length(v1) > 1) 
        stop("v1 has to be a scalar")
    
    if (length(theta) > 1) 
        stop("theta has to be a scalar")
    
    
    print("Loading Data")
    #################################################################
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    
    qq <- ncol(X)
    nn <- nrow(X)
    
    if(nn != dim(Y)[1])
        stop("length of Y for each voxel has to be nrow(X)")
    
    n <- nn / 2
    q <- qq / 2
    x.cv <- X[1:(n), 1:q]
    N <- dim(Y)[2] # matrix size
    # XtX <- crossprod(X)
    XtX.cv <- crossprod(x.cv)
    Yre <- Y[1:(n), , ]
    Yim <- Y[-c(1:(n)), , ]
    L <- length(v0s)
    omega1 <- 2 * v1
    lam1 <- 2i * v1 * slab.cor
    b1 <- Mod(lam1) ^ 2 / (Mod(omega1) ^ 2 - Mod(lam1) ^ 2)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sigmas <- rep(0, L)
    betas <- array(0, c(L, qq, N, N))  # gamma coef in the paper
    beta_k <- array(0, c(qq, N, N))
    beta_new <- array(0, c(qq, N, N))
    postprobs <- array(0, c(L, q, N, N)) 
    p_star_mat <- array(0, c(q, N, N))
    b_star <- array(0, c(q, N, N))
    # Bmat <- array(0, c(qq, qq, N, N))
    log_gfcn <- rep(0, L)
    counts <- rep(0, L)
    thetas <- rep(0, L)
    
    
    print("General EMcplx spikecir algo Begins")
    #-------------------------------
    # EM algorithm
    #------------------------------- 
    for (m in 1:L) {
        v0 <- v0s[m]
        print(paste("v0 = ", v0))
        omega0 <- 2 * v0
        lam0 <- 2i * v0 * spike.cor
        b0 <- Mod(lam0) ^ 2 / (Mod(omega0) ^ 2 - Mod(lam0) ^ 2)
        eps <- epsilon + 1
        count <- 1
        
        # Initialization of parameters
        #----------------------------
        theta_k <- 0.5
        sigma_k <- 1
        sigma_new <- sigma_k
        theta_new <- theta_k
        
        #-------------------------------
        # Starting values
        #-------------------------------
        Beta_init <- array(0, c(qq, N, N))
        for (j in 1:N) {
            for (i in 1:N) {
                # initial values: Start with LSE
                lmout <- lm(Y[, i, j] ~ 0 + X)
                Beta_init[, i, j] <- lmout$coeff
                beta_k[, i, j] <- Beta_init[, i, j]
                beta_new[, i, j] <- beta_k[, i, j]
            }
        }
        while(eps > epsilon) {
            if (count == 2000) break
            ss_lik <- 0
            ss_D <- 0
            ss_G <- 0
            for (j in 1:N) {
                for (i in 1:N) {
                    # ******** E-step ************ #
                    # Update inclusion probabilities
                        
                    p_star_v <- E_beta_binom_spikecir_cpp(beta_k[, i, j], sigma_k, 
                                                      theta_k, as.matrix(omega0), 
                                                      as.matrix(omega1), 
                                                      as.matrix(lam0), 
                                                      as.matrix(lam1))
                    p_star_mat[, i, j] <- p_star_v
                    # print(p_star_mat[, i, j])
                    
                    # ******** M-step for beta *** #
                    XtYre <- crossprod(x.cv, Yre[, i, j])
                    XtYim <- crossprod(x.cv, Yim[, i, j])
                    Dv <- create_Dmat_noncir_cpp(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    # print(Dv)
                    Gv <- create_Gmat_noncir_cpp(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv_re <- Re(Gv)
                    Gv_im <- Im(Gv)
                    # print(Gv)
                    # complex-valued reg coef
                    ga_v <- M_beta_noncir_cpp(XtYre, XtYim, XtX.cv, Dv, Gv_re, Gv_im)
                    # print(ga_v)
                    # ga_v <- ga_v$beta_cx
                    # save re and im parts separately in real representation
                    beta_k[, i, j] <- c(Re(ga_v), Im(ga_v))
                    
                    # compute sum of squars for sigma 
                    Y.cplx <- Yre[, i, j] + Yim[, i, j] * 1i
                    # print(head(Y.cplx, 2))
                    res.cv <- Y.cplx - x.cv %*% ga_v
                    # print(head(x.cv %*% ga_v, 2))
                    # print(head(res.cv, 2))
                    ss_lik <- ss_lik + Re(t(Conj(res.cv)) %*% res.cv)
                    ss_D <- ss_D + Re(t(Conj(ga_v)) %*% Dv %*% ga_v)
                    ss_G <- ss_G + Re(t(ga_v) %*% Gv %*% ga_v)
                }
            }
            # ******** M-step for theta and sigma ************ #
            theta_k <- M_theta_cpp(p_star_mat, a_th, b_th)
            sigma_k <- M_sigma_sig_noncir_cpp(ss_lik, ss_D, ss_G, a_sig=1/2, b_sig=1/2, N, n, q) #check
            eps <- sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            
            beta_new <- beta_k
            sigma_new <- sigma_k
            theta_new <- theta_k
            count <- count + 1
            print(paste("Epsilon", eps, " count", count))
            #             print(paste("sigma_k", sigma_k))
            #             print(paste("theta_k", theta_k))
        }
        
        #  Save values 
        # -----------------------------------------------------
        postprobs[m, , , ] <- p_star_mat
        betas[m, , , ] <- beta_new
        sigmas[m] <- sigma_new
        thetas[m] <- theta_new
        counts[m] <- count
        Active <- (p_star_mat > thre)
        #         log_g_sig_noncir <- function(Active, X, Y, a_sig, b_sig,
        #                                      lam0, lam1, omega0, omega1, b0, b1, a, b)
        log_gfcn[m] <- log_g_sig_noncir(Active, X, Y, a_sig=1/2, b_sig=1/2, p_star_mat,
                                        lam0, lam1, omega0, omega1, b0, b1, a_th, b_th)
        # Ymat <- matrix(Y, 400, 2304)
        # log_gfcn[m] <- log_g_sig_noncir_cpp(Active, X, Ymat, 1/2, 1/2, p_star_mat,
        #                                 lam0, lam1, omega0, omega1, b0, b1, a_th, b_th)
    }
    print("Done!")
    return(list(pprob = postprobs, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts, Beta_init = Beta_init))
}

EM_spikecir_cpp_full <- function(X, Y, v0s, v1, epsilon, theta, a_th, b_th, thre, 
                            spike.cor, slab.cor) {
    # X: n x qq matrix 
    # Y: n x N x N Array
    # v0: vector of length L
    # v1: scalar
    # beta_init: qq x 1 vector
    # sigma_init: scalar
    # epsilon: scalar
    # theta: scalar
    # a: scalar
    # b: scalar
    # va_g: scalar
    #################################################################
    if (missing(theta)){theta = 0.5}
    
    if (class(X) != "matrix") {
        tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("X must be a matrix or able to be coerced to a matrix")
    }
    
    if (class(Y) != "array") {
        tmp <- try(Y <- as.array(Y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("Y must be an array or able to be coerced to an array")
    }
    
    if (sum(v1<v0s) > 0)
        stop("v1 has to be larger than v0")
    
    if (length(a_th) > 1) 
        stop("a has to be a scalar")
    
    if (length(b_th) > 1) 
        stop("b has to be a scalar")
    
    if (length(v1) > 1) 
        stop("v1 has to be a scalar")
    
    if (length(theta) > 1) 
        stop("theta has to be a scalar")

    #################################################################
    Ymat <- matrix(Y, dim(Y)[1], dim(Y)[2] ^ 2)
    EM_spikecir_main_cpp(X, Y, Ymat, v0s, v1, epsilon, theta, a_th, b_th, 
                         thre, spike.cor, slab.cor, log_g_sig_noncir)
}











