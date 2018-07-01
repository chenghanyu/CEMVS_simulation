# M-step for the standard deviation (sigma) for EMVSmain_sig


M_sigma_sig <- function(Y, X, beta_k, d_star, eta, lambda) {
    # Y: 3-dim array
    # beta_k: 3-dim array
    # d_star: 3-dim array
    NN = dim(Y)[2] ^ 2
    TT = nrow(X)
    q = ncol(X)
    den = NN * (TT + q) + 2 + eta
    num = 0
    for (j in 1:N) {
        for (i in 1:N) {
            num = num + crossprod(Y[, i, j] - X %*% beta_k[, i, j]) + 
                crossprod(sqrt(diag(as.numeric(d_star[, i, j]), ncol = q, 
                                    nrow = q)) %*% beta_k[, i, j])
        }
    }

    sig <- as.numeric(sqrt(num / den))
    return(sig)
}

M_sigma_sig_cp <- function(Y, X, beta_k, d_star, eta, lambda) {
    N = dim(Y)[2]
    TT = nrow(X)
    qq = ncol(X)
    den = N ^ 2 * (TT + qq) + 2 + eta
    num = 0
    for (j in 1:N) {
        for (i in 1:N) {
            num = num + crossprod(Y[, i, j] - X %*% beta_k[, i, j]) + 
                crossprod(sqrt(diag(as.numeric(d_star[, i, j]), ncol = qq, 
                                    nrow = qq)) %*% beta_k[, i, j])
        }
    }
    sig <- as.numeric(sqrt(num / den))
    return(sig)
}