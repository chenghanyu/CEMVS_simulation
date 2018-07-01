# log g function for EMVSmain_sig

# source("./M_beta_noncir.R")
# source("./log_prior.R")

log_g_sig_cp_noncir <- function (Gamma, X, Y, eta, lambda, D, v0, v1, a, b) {
    # D is a 2p by 2p matrix
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    pp <- ncol(X)
    p <- pp / 2
    n <- nrow(X)
    N <- dim(Y)[2]
    Ssq_sum <- 0
    XX_sum <- 0
    dd_sum <- 0
    if (v0 > 0) {
        if (pp > 0) {
            X_gamma = data.matrix(X)
            XtX = t(X_gamma) %*% X_gamma
            for (j in 1:N) {
                for (i in 1:N) {
                    XtY = t(X_gamma) %*% Y[ , i, j]
                    dd <- (1 - Gamma[, i, j]) * v0 + Gamma[, i, j] * v1
                    X_tilde = rbind(X_gamma, sqrt(1 / dd))
                    aux = M_beta(XtY, X_gamma, Y[, i, j], XtX, 
                                 d_star = rep(1 / dd, 2))
                    Ssq = t(Y[, i, j]) %*% Y[, i, j] - t(XtY) %*% aux
                    Ssq_sum = Ssq_sum + Ssq
                    xx = (-0.5) * as.numeric(
                        determinant(t(as.matrix(X_tilde)) %*% as.matrix(X_tilde), 
                                    logarithm = TRUE)$modulus)
                    XX_sum = XX_sum + xx
                    ddd =  -(1/2) * sum(log(rep(dd, 2)))
                    dd_sum = dd_sum + ddd
                }
            }
            value = XX_sum + dd_sum  - (0.5) * (N ^ 2 * n + eta) * 
                log(eta * lambda + Ssq_sum)
            value = value + log_prior(Gamma, a, b, length(Gamma))
        }
    }
    if (v0 == 0) {
        for (j in 1:N) {
            for (i in 1:N) {
                p <- sum(Gamma[, i, j])
                if (p > 0) {
                    Dinv <- chol2inv(chol(D))
                    X_gamma <- data.matrix(X[, rep(Gamma[, i, j], 2) == 1])
                    XtY <- crossprod(X_gamma, Y[, i, j])
                    XtX <- crossprod(X_gamma)
                    aux <- M_beta_cp_noncir(XtY, XtX, Dinv / v1)
                    Ssq <- crossprod(Y[, i, j])- crossprod(XtY, aux)
                    Ssq_sum <- Ssq_sum + Ssq
                    det.value <- (-0.5) * as.numeric(
                        determinant(XtX + (Dinv / v1), logarithm = TRUE)$modulus)
                    XX_sum <- XX_sum + det.value
                    ddd <- (-0.5) * as.numeric(determinant(v1 * D, logarithm = TRUE)$modulus)
                    dd_sum <- dd_sum + ddd
                } else {
                    Ssq <- crossprod(Y[, i, j])
                    Ssq_sum <- Ssq_sum + Ssq
                }
            }
        }
        value <- XX_sum + dd_sum - (0.5) * (N ^ 2 * n + eta) * 
            log(eta * lambda + Ssq_sum)
        value <- value + log_prior(Gamma, a, b, length(Gamma))
    }
    return(value)
}

log_g_sig_noncir <- function(Active, X, Y, a_sig, b_sig, p_star_mat,
                             lam0, lam1, omega0, omega1, b0, b1, a, b, v0) {
    # D is a 2p by 2p matrix
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    qq <- ncol(X)
    q <- qq / 2
    nn <- nrow(X)
    n <- nn / 2
    N <- dim(Y)[2]
    Ssq_sum <- 0
    XX_sum <- 0
    dd_sum <- 0
    Yre <- Y[1:(n), , ]
    Yim <- Y[-c(1:(n)), , ]
    if (v0 == 0) {
        b1 <- Mod(lam1) ^ 2 / (Mod(omega1) ^ 2 - Mod(lam1) ^ 2)
        if (length(dim(Active)) == 2) Active = array(Active, n, N, N)
        for (j in 1:N) {
            for (i in 1:N) {
                no.active <- sum(Active[, i, j])
                if (no.active > 0) {
                    x <- X[1:n, 1:q]
                    x_act <- data.matrix(as.matrix(x)[, Active[, i, j] == 1])
                    X_act <- data.matrix(X[, rep(Active[, i, j], 2) == 1])

                    XtYre <- crossprod(x_act, Yre[, i, j])
                    XtYim <- crossprod(x_act, Yim[, i, j])
                    Dv <- create_Dmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star = p_star_mat[,i, j])
                    Gv <- create_Gmat_noncir(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star = p_star_mat[,i, j])

                    XtX.cv <- crossprod(x_act)
                    regcoef <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Re(Gv), Im(Gv))
                    XtY <- crossprod(X_act, Y[, i, j])
                    Ssq <- crossprod(Y[, i, j]) - crossprod(XtY, c(Re(regcoef), Im(regcoef)))
                    Ssq_sum <- Ssq_sum + Ssq
                    XtX <- crossprod(X_act)
                    Sig.real <- cplx_to_real_cov(Omega = omega1, Lam = lam1)
                    det.value <- (-0.5) * as.numeric(
                        determinant(XtX + chol2inv(chol(Sig.real)), 
                                    logarithm = TRUE)$modulus)
                    print(paste("det.value  ", det.value))
                    XX_sum <- XX_sum + det.value
                    ddd <- (-0.5) * as.numeric(determinant(Sig.real, 
                                                           logarithm = TRUE)$modulus)
                    dd_sum <- dd_sum + ddd
                } else {
                    Ssq <- crossprod(Y[, i, j])
                    Ssq_sum <- Ssq_sum + Ssq
                }
            }
        }
        print(paste("Ssq_sum  ",Ssq_sum))
        print(paste("XX_sum  ", XX_sum))
        print(paste("dd_sum  ", dd_sum))
        value <- XX_sum + dd_sum - ((0.5) * (N ^ 2 * nn) + a_sig) * 
            log(b_sig + Ssq_sum)
        value <- value + log_prior(Active, a, b, length(Active))
        print(paste("log_prior  ", log_prior(Active, a, b, length(Active))))
    }
    return(value)
}

log_g_sig_noncir <- function(Active, X, Ycube, a_sig, b_sig, p_star_mat,
                             lam0, lam1, omega0, omega1, b0, b1, a, b, v0) {
    # D is a 2p by 2p matrix
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    qq <- ncol(X)
    q <- qq / 2
    nn <- nrow(X)
    n <- nn / 2
    N <- dim(Ycube)[2]
    Ssq_sum <- 0
    XX_sum <- 0
    dd_sum <- 0
    Yre <- Ycube[1:(n), , ]
    Yim <- Ycube[-c(1:(n)), , ]
    if (v0 == 0) {
        b1 <- Mod(lam1) ^ 2 / (Mod(omega1) ^ 2 - Mod(lam1) ^ 2)
        if (length(dim(Active)) == 2) Active <- array(Active, dim = c(q, N, N))
        if (length(dim(p_star_mat)) == 2) p_star_mat <- array(p_star_mat, dim = c(q, N, N))
        for (j in 1:N) {
            for (i in 1:N) {
                no.active <- sum(Active[, i, j])
                if (no.active > 0) {
                    x <- X[1:n, 1:q]
                    x_act <- data.matrix(as.matrix(x)[, Active[, i, j] == 1])
                    X_act <- data.matrix(X[, rep(Active[, i, j], 2) == 1])

                    XtYre <- crossprod(x_act, Yre[, i, j])
                    XtYim <- crossprod(x_act, Yim[, i, j])
                    Dv <- create_Dmat_noncir_cpp(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_mat[,i, j])
                    Gv <- create_Gmat_noncir_cpp(lam0, lam1, omega0, omega1, b0, b1, 
                                             p_star_mat[,i, j])
                    XtX.cv <- crossprod(x_act)

                    regcoef <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Re(Gv), Im(Gv))
                    XtY <- crossprod(X_act, Ycube[, i, j])
                    Ssq <- crossprod(Ycube[, i, j]) - crossprod(XtY, c(Re(regcoef), Im(regcoef)))
                    Ssq_sum <- Ssq_sum + Ssq
                    XtX <- crossprod(X_act)
                    Sig.real <- cplx_to_real_cov(Omega = omega1, Lam = lam1)
                    det.value <- (-0.5) * as.numeric(
                        determinant(XtX + chol2inv(chol(Sig.real)), 
                                    logarithm = TRUE)$modulus)
                    XX_sum <- XX_sum + det.value
                    ddd <- (-0.5) * as.numeric(determinant(Sig.real, 
                                                           logarithm = TRUE)$modulus)
                    dd_sum <- dd_sum + ddd
                } else {
                    Ssq <- crossprod(Ycube[, i, j])
                    Ssq_sum <- Ssq_sum + Ssq
                }
            }
        }
        value <- XX_sum + dd_sum - ((0.5) * (N ^ 2 * nn) + a_sig) * 
            log(b_sig + Ssq_sum)
        value <- value + log_prior(Active, a, b, length(Active))
    }
    return(value)
}

log_g_sig_noncir_c <- function(Active, X, Y, a_sig, b_sig, p_star_mat,
                             lam0, lam1.mat, omega0, omega1, b0, b1, a, b, v0=0) {  
    # D is a 2p by 2p matrix
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    qq <- ncol(X)
    q <- qq / 2
    nn <- nrow(X)
    n <- nn / 2
    N <- dim(Y)[2]
    Ssq_sum <- 0
    XX_sum <- 0
    dd_sum <- 0
    Yre <- Y[1:(n), , ]
    Yim <- Y[-c(1:(n)), , ]
    if (v0 == 0) {

        for (j in 1:N) {
            for (i in 1:N) {
                n.active <- sum(Active[, i, j])
                if (n.active > 0) {
                    x <- X[1:n, 1:q]
                    x_act <- data.matrix(as.matrix(x)[, Active[, i, j] == 1])
                    X_act <- data.matrix(X[, rep(Active[, i, j], 2) == 1])
                    XtYre <- crossprod(x_act, Yre[, i, j])
                    XtYim <- crossprod(x_act, Yim[, i, j])
                    b1 <- Mod(lam1.mat[i, j]) ^ 2 / (Mod(omega1) ^ 2 - Mod(lam1.mat[i, j]) ^ 2)
                    Dv <- create_Dmat_noncir(lam0, lam1.mat[i, j], omega0, omega1, b0, b1, 
                                             p_star = p_star_mat[,i, j])
                    Gv <- create_Gmat_noncir(lam0, lam1.mat[i, j], omega0, omega1, b0, b1, 
                                             p_star = p_star_mat[,i, j])
                    XtX.cv <- crossprod(x_act)
                    regcoef <- M_beta_noncir(XtYre, XtYim, XtX.cv, Dv, Re(Gv), Im(Gv))
                    XtY <- crossprod(X_act, Y[, i, j])
                    Ssq <- crossprod(Y[, i, j])- crossprod(XtY, c(Re(regcoef), Im(regcoef)))
                    Ssq_sum <- Ssq_sum + Ssq
                    
                    XtX <- crossprod(X_act)
                    Sig.real <- cplx_to_real_cov(Omega = omega1, Lam = lam1.mat[i, j])
                    det.value <- (-0.5) * as.numeric(
                        determinant(XtX + chol2inv(chol(Sig.real)), 
                                    logarithm = TRUE)$modulus)
                    XX_sum <- XX_sum + det.value
                    ddd <- (-0.5) * as.numeric(determinant(Sig.real, 
                                                           logarithm = TRUE)$modulus)
                    dd_sum <- dd_sum + ddd
                } else {
                    Ssq <- crossprod(Y[, i, j])
                    Ssq_sum <- Ssq_sum + Ssq
                }
            }
        }
        value <- XX_sum + dd_sum - ((0.5) * (N ^ 2 * nn) + a_sig) * 
                    log(b_sig + Ssq_sum)
        value <- value + log_prior(Active, a, b, length(Active))
    }
    return(value)
}
