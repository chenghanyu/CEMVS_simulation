# EMVS noncircular one theta one sig main function 
# This is the most General version including matrix G

cplx_to_real_cov <- function(Omega, Lam) {
    Vrr <- (0.5) * Re(Omega + Lam)
    Vii <- (0.5) * Re(Omega - Lam)
    Vri <- (0.5) * Im(-Omega + Lam)
    Vir <- (0.5) * Im(Omega + Lam)
    return(rbind(cbind(Vrr, Vri), cbind(Vir, Vii)))
}

real_to_cplx_cov <- function(Vrr, Vri, Vir, Vii) {
    Omega <- Vrr + Vii + (Vir - Vri) * 1i
    Lam <- Vrr - Vii + (Vir + Vri) * 1i
    return(list(Omega = Omega, Lam = Lam))
}


EM_noncir <- function(X, Y, v0s, v1, epsilon, theta, a_th, b_th, v1_g, thre, 
                                cor.ri) {
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
    if (missing(v1_g)){v1_g=v1}
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
    
    if (length(v1_g) >1) 
        stop("v1_g has to be a scalar")

    
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
    XtX.cv <- crossprod(x.cv)
    Yre <- Y[1:(n), , ]
    Yim <- Y[-c(1:(n)), , ]
    L <- length(v0s)
    omega1 <- 2 * v1
    lam1 <- 2i * v1 * cor.ri
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
    


    print("General EMcplx noncir algo Begins")
    #-------------------------------
    # EM algorithm
    #------------------------------- 
    for (m in 1:L) {
        v0 <- v0s[m]
        print(paste("v0 = ", v0))
        omega0 <- 2 * v0
        lam0 <- 2i * v0 * cor.ri
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
                    p_star_v <- E_beta_binom_noncir(beta_k[, i, j], 
                                                     sigma_k, theta_k, 
                                                     omega0, omega1, 
                                                     lam0, lam1)
                    p_star_mat[, i, j] <- p_star_v

                    # ******** M-step for beta *** #
                    XtYre <- crossprod(x.cv, Yre[, i, j])
                    XtYim <- crossprod(x.cv, Yim[, i, j])
                    Dv <- create_Dmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv <- create_Gmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv_re <- Re(Gv)
                    Gv_im <- Im(Gv)
                    # complex-valued reg coef
                    ga_v <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Gv_re, Gv_im)
                    
                    # save re and im parts separately in real representation
                    beta_k[, i, j] <- c(Re(ga_v), Im(ga_v))
                    
                    # compute sum of squars for sigma 
                    Y.cplx <- Yre[, i, j] + Yim[, i, j] * 1i
                    res.cv <- Y.cplx - x.cv %*% ga_v
                    ss_lik <- ss_lik + Re(t(Conj(res.cv)) %*% res.cv)
                    ss_D <- ss_D + Re(t(Conj(ga_v)) %*% Dv %*% ga_v)
                    ss_G <- ss_G + Re(t(ga_v) %*% Gv %*% ga_v)
                }
            }
            # ******** M-step for theta and sigma ************ #
            theta_k <- M_theta_cp(p_star_mat, a_th, b_th)
            sigma_k <- M_sigma_sig_noncir(ss_lik, ss_D, ss_G, a_sig=1/2, b_sig=1/2, N, n, q) #check
            eps <- sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            
            beta_new <- beta_k
            sigma_new <- sigma_k
            theta_new <- theta_k
            count <- count + 1
            print(paste("Epsilon", eps, " count", count))
        }
        
        #  Save values 
        # -----------------------------------------------------
        postprobs[m, , , ] <- p_star_mat
        betas[m, , , ] <- beta_new
        sigmas[m] <- sigma_new
        thetas[m] <- theta_new
        counts[m] <- count
        Active <- (p_star_mat > thre)
        log_gfcn[m] <- log_g_sig_noncir(Active, X, Y, a_sig=1/2, b_sig=1/2, p_star_mat,
                                        lam0, lam1, omega0, omega1, b0, b1, a_th, b_th)
    }
    print("Done!")
    return(list(pprob = postprobs, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts, Beta_init = Beta_init))
}

EM_spikecir <- function(X, Y, v0s, v1, epsilon, theta, a_th, b_th, v1_g, thre, 
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
    if (missing(v1_g)){v1_g=v1}
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
    
    if (length(v1_g) >1) 
        stop("v1_g has to be a scalar")
    
    
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
                    p_star_v <- E_beta_binom_spikecir(beta_k[, i, j], sigma_k, 
                                                      theta_k, omega0, omega1, 
                                                      lam0, lam1)
                    p_star_mat[, i, j] <- p_star_v
                    
                    # ******** M-step for beta *** #
                    XtYre <- crossprod(x.cv, Yre[, i, j])
                    XtYim <- crossprod(x.cv, Yim[, i, j])
                    Dv <- create_Dmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv <- create_Gmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv_re <- Re(Gv)
                    Gv_im <- Im(Gv)
                    # complex-valued reg coef
                    ga_v <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Gv_re, Gv_im)
                    
                    # save re and im parts separately in real representation
                    beta_k[, i, j] <- c(Re(ga_v), Im(ga_v))
                    
                    # compute sum of squars for sigma 
                    Y.cplx <- Yre[, i, j] + Yim[, i, j] * 1i
                    res.cv <- Y.cplx - x.cv %*% ga_v
                    ss_lik <- ss_lik + Re(t(Conj(res.cv)) %*% res.cv)
                    ss_D <- ss_D + Re(t(Conj(ga_v)) %*% Dv %*% ga_v)
                    ss_G <- ss_G + Re(t(ga_v) %*% Gv %*% ga_v)
                }
            }
            # ******** M-step for theta and sigma ************ #
            theta_k <- M_theta_cp(p_star_mat, a_th, b_th)
            sigma_k <- M_sigma_sig_noncir(ss_lik, ss_D, ss_G, a_sig=1/2, b_sig=1/2, N, n, q) #check
            eps <- sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            
            beta_new <- beta_k
            sigma_new <- sigma_k
            theta_new <- theta_k
            count <- count + 1
            print(paste("Epsilon", eps, " count", count))
        }
        
        #  Save values 
        # -----------------------------------------------------
        postprobs[m, , , ] <- p_star_mat
        betas[m, , , ] <- beta_new
        sigmas[m] <- sigma_new
        thetas[m] <- theta_new
        counts[m] <- count
        Active <- (p_star_mat > thre)
        log_gfcn[m] <- log_g_sig_noncir(Active, X, Y, a_sig=1/2, b_sig=1/2, p_star_mat,
                                        lam0, lam1, omega0, omega1, b0, b1, a_th, b_th)
    }
    print("Done!")
    return(list(pprob = postprobs, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts, Beta_init = Beta_init))
}

EM_spikecir_c <- function(X, Y, v0s, v1, epsilon, theta, a_th, b_th, v1_g, thre, 
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
    if (missing(v1_g)){v1_g=v1}
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
    
    if (length(v1_g) >1) 
        stop("v1_g has to be a scalar")
    
    
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
        lam1.mat <- matrix(0, N, N)
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
                    p_star_v <- E_beta_binom_spikecir(beta_k[, i, j], sigma_k, 
                                                      theta_k, omega0, omega1, 
                                                      lam0, lam1)
                    p_star_mat[, i, j] <- p_star_v
                    
                    lam1 <- 2i * v1 * slab.cor * p_star_v
                    lam1.mat[i, j] <- lam1
                    
                    # ******** M-step for beta *** #
                    XtYre <- crossprod(x.cv, Yre[, i, j])
                    XtYim <- crossprod(x.cv, Yim[, i, j])
                    Dv <- create_Dmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv <- create_Gmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_v)
                    Gv_re <- Re(Gv)
                    Gv_im <- Im(Gv)
                    # complex-valued reg coef
                    ga_v <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Gv_re, Gv_im)
                    
                    # save re and im parts separately in real representation
                    beta_k[, i, j] <- c(Re(ga_v), Im(ga_v))
                    
                    # compute sum of squars for sigma 
                    Y.cplx <- Yre[, i, j] + Yim[, i, j] * 1i
                    res.cv <- Y.cplx - x.cv %*% ga_v
                    ss_lik <- ss_lik + Re(t(Conj(res.cv)) %*% res.cv)
                    ss_D <- ss_D + Re(t(Conj(ga_v)) %*% Dv %*% ga_v)
                    ss_G <- ss_G + Re(t(ga_v) %*% Gv %*% ga_v)
                }
            }
            # ******** M-step for theta and sigma ************ #
            theta_k <- M_theta_cp(p_star_mat, a_th, b_th)
            sigma_k <- M_sigma_sig_noncir(ss_lik, ss_D, ss_G, a_sig=1/2, b_sig=1/2, N, n, q) #check
            eps <- sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            
            beta_new <- beta_k
            sigma_new <- sigma_k
            theta_new <- theta_k
            count <- count + 1
            print(paste("Epsilon", eps, " count", count))
        }
        
        #  Save values 
        # -----------------------------------------------------
        postprobs[m, , , ] <- p_star_mat
        betas[m, , , ] <- beta_new
        sigmas[m] <- sigma_new
        thetas[m] <- theta_new
        counts[m] <- count
        Active <- (p_star_mat > thre)
        log_gfcn[m] <- log_g_sig_noncir_c(Active, X, Y, a_sig=1/2, b_sig=1/2, p_star_mat,
                                        lam0, lam1.mat, omega0, omega1, b0, b1, a_th, b_th)
    }
    print("Done!")
    return(list(pprob = postprobs, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts, Beta_init = Beta_init))
}


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


















