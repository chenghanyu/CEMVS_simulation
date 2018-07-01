# M-step for the standard deviation (sigma) for EMVSmain_sig

M_sigma_sig_cp_noncir <- function(Y, X, beta_k, B, eta, lambda) {
    # B: an 4-dim array (for each voxel there is a 2p by 2p matrix)
    N <- dim(Y)[2]
    TT <- nrow(X)
    qq <- ncol(X)
    den <- N ^ 2 * (TT + qq) + 2 + eta
    num <- 0
    for (j in 1:N) {
        for (i in 1:N) {
            num <- num + crossprod(Y[, i, j] - X %*% beta_k[, i, j]) + 
                crossprod(chol(B[ , , i, j]) %*% beta_k[, i, j])
        }
    }
    sig <- as.numeric(sqrt(num / den))
    return(sig)
}

M_sigma_sig_noncir <- function(ss_lik, ss_D, ss_G, a_sig, b_sig, N, n, q) {
    den <- N ^ 2 * (n + q) + 1 + a_sig
    num <- ss_lik/2 + ss_D - ss_G + b_sig
    sig <- as.numeric(sqrt(num / den))
    return(sig)
}


# same sigma update
M_sigma_sig_cp_spikecir <- function(Y, X, beta_k, B, eta, lambda) {
    # B: an 4-dim array (for each voxel there is a 2p by 2p matrix)
    N <- dim(Y)[2]
    TT <- nrow(X)
    qq <- ncol(X)
    den <- N ^ 2 * (TT + qq) + 2 + eta
    num <- 0
    for (j in 1:N) {
        for (i in 1:N) {
            num <- num + crossprod(Y[, i, j] - X %*% beta_k[, i, j]) + 
                crossprod(chol(B[ , , i, j]) %*% beta_k[, i, j])
        }
    }
    sig <- as.numeric(sqrt(num / den))
    return(sig)
}






