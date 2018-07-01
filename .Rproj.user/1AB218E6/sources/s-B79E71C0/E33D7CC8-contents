# M-step for regression coefficients 

M_beta <- function(XtY, X, Y, XtX, d_star) {
    n = nrow(X)
    q = ncol(X)
    if (q > n) {
        X_star = t(t(X) / (tausq <- as.numeric(d_star)))
        A = solve(diag(n) + X_star %*% t(X), tol = 1e-400)
        result = ((1 / tausq) * diag(length(tausq)) - t(X_star) %*% A %*% X_star)
        result = result %*% XtY
        return(result)
    } else {
        Psi = XtX
        diag(Psi) = diag(Psi) + as.numeric(d_star)
        Psi = chol2inv(chol(Psi))
        Psi = Psi %*% XtY
        return(Psi)
    }
}

M_beta_cp <- function(XtY, X, Y, XtX, d_star) {
    # d_star: dim = q = qq / 2
    n = nrow(X)
    qq = ncol(X)
    if (qq > n) {
        X_star = t(t(X) / (tausq <- rep(as.numeric(d_star), 2)))
        A = solve(diag(n) + X_star %*% t(X), tol = 1e-400)
        result = ((1 / tausq) * diag(length(tausq)) - t(X_star) %*% A %*% X_star)
        result = result %*% XtY
        return(result)
    } else {
        Psi = XtX
        diag(Psi) = diag(Psi) + rep(as.numeric(d_star), 2)
        Psi = chol2inv(chol(Psi))
        Psi = Psi %*% XtY
        return(Psi)
    }
}
