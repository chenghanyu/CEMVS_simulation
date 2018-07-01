# log g function for EMVSmain_sig

log_g_sig <- function (Gamma, X, Y, eta, lambda, v0, v1, a, b) {
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    p = ncol(X)
    n = nrow(X)
    N = dim(Y)[2]
    Ssq_sum = 0
    XX_sum = 0
    dd_sum = 0
    if (v0 > 0) {
        if (p > 0) {
            X_gamma = data.matrix(X)
            XtX = t(X_gamma) %*% X_gamma
            # X_gamma = data.matrix(X[, gamma == 1])
            for (j in 1:N) {
                for (i in 1:N) {
                    XtY = t(X_gamma) %*% Y[ , i, j]
                    dd <- (1 - Gamma[, i, j]) * v0 + Gamma[, i, j] * v1
                    X_tilde = rbind(X_gamma, sqrt(1 / dd))
                    aux = M_beta(XtY, X_gamma, Y[, i, j], XtX, d_star = 1 / dd)
                    Ssq = t(Y[, i, j]) %*% Y[, i, j] - t(XtY) %*% aux
                    Ssq_sum = Ssq_sum + Ssq
                    xx = (-0.5) * as.numeric(
                        determinant(t(as.matrix(X_tilde)) %*% as.matrix(X_tilde), 
                                    logarithm = TRUE)$modulus)
                    XX_sum = XX_sum + xx
                    ddd =  -(1/2) * sum(log(dd))
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
                p = sum(Gamma[, i, j])
                if (p > 0) {
                    X_gamma = data.matrix(X[, Gamma[, i, j] == 1])
                    X_tilde = rbind(X_gamma, sqrt(1 / v1) * diag(p))
                    XtY = crossprod(X_gamma, Y[, i, j])
                    XtX = crossprod(X_gamma)
                    aux = M_beta(XtY, X_gamma, Y[, i, j], XtX, rep(1 / v1, p))
                    Ssq = crossprod(Y[, i, j])- t(XtY) %*% aux
                    Ssq_sum = Ssq_sum + Ssq
                    xx = (-0.5) * as.numeric(
                        determinant(t(as.matrix(X_tilde)) %*% as.matrix(X_tilde), 
                                    logarithm = TRUE)$modulus)
                    XX_sum = XX_sum + xx
                    ddd =  -p/2 * log(v1)
                    dd_sum = dd_sum + ddd
                } else {
                    Ssq <- crossprod(Y[, i, j])
                    Ssq_sum = Ssq_sum + Ssq
                }
            }
        }
        value = XX_sum + dd_sum - (0.5) * (N ^ 2 * n + eta) * 
            log(eta * lambda + Ssq_sum)
        value = value + log_prior(Gamma, a, b, length(Gamma))
    }
    return(value)
}

log_g_sig_cp <- function (Gamma, X, Y, eta, lambda, v0, v1, a, b) {
    # Assume R = I
    # Gamma: q * N * N
    # Y: 3d array
    pp = ncol(X)
    p = pp / 2
    n = nrow(X)
    N = dim(Y)[2]
    Ssq_sum = 0
    XX_sum = 0
    dd_sum = 0
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
        # index = p_star[, i, j] > thre
        # gamma = as.numeric(index)
        for (j in 1:N) {
            for (i in 1:N) {
                p = sum(Gamma[, i, j])
                if (p > 0) {
                    X_gamma = data.matrix(X[, rep(Gamma[, i, j], 2) == 1])
                    X_tilde = rbind(X_gamma, sqrt(1 / v1) * diag(2*p))
                    XtY = t(X_gamma) %*% Y[, i, j]
                    XtX = t(X_gamma) %*% X_gamma
                    aux = M_beta(XtY, X_gamma, Y[, i, j], XtX, rep(1 / v1, 2*p))
                    Ssq = t(Y[, i, j]) %*% Y[, i, j] - t(XtY) %*% aux
                    Ssq_sum = Ssq_sum + Ssq
                    xx = (-0.5) * as.numeric(
                        determinant(t(as.matrix(X_tilde)) %*% as.matrix(X_tilde), 
                                    logarithm = TRUE)$modulus)
                    XX_sum = XX_sum + xx
                    ddd =  -p/2 * log(v1)
                    dd_sum = dd_sum + ddd
                } else {
                    Ssq <- t(Y[, i, j]) %*% Y[, i, j]
                    Ssq_sum = Ssq_sum + Ssq
                }
            }
        }
        value = XX_sum + dd_sum - (0.5) * (N ^ 2 * n + eta) * 
            log(eta * lambda + Ssq_sum)
        value = value + log_prior(Gamma, a, b, length(Gamma))
    }
    return(value)
}
