# M-step for regression coefficients (beta, gamma)

M_beta_cp_noncir <- function(XtY, XtX, Bv) {
    # d_star: dim = q = qq / 2
    # n <- nrow(X)
    # qq <- ncol(X)
    if (ncol(X) > nrow(X)) {
        X_star = t(t(X) / (tausq <- rep(as.numeric(d_star), 2)))
        A = solve(diag(n) + X_star %*% t(X), tol = 1e-400)
        result = ((1 / tausq) * diag(length(tausq)) - t(X_star) %*% A %*% X_star)
        result = result %*% XtY
        return(result)
    } else {
        Psi <- chol2inv(chol(XtX + Bv))
        return(Psi %*% XtY)
    }
}


create_Dmat_noncir <- function(lam0, lam1, omega0, omega1, b0, b1, p_star) {
    # p_star: vector fof length p
    # lam0, lam1, omega0, omega1: element in Omega and Lambda matrix for 
    #                             spike and slab parts
    D <- (1 / omega0) * (1 + b0) * (1 - p_star) + 
        (1 / omega1) * (1 + b1) * p_star
    return(diag(D, length(p_star)))
    
}

create_Gmat_noncir <- function(lam0, lam1, omega0, omega1, b0, b1, p_star) {
    # p_star: vector fof length p
    # lam0, lam1, omega0, omega1: element in Omega and Lambda matrix for 
    #                             spike and slab parts
    
    #     b0 <- Mod(lam0) ^ 2 / (Mod(omega0) ^ 2 - Mod(lam0) ^ 2)
    #     b1 <- Mod(lam1) ^ 2 / (Mod(omega1) ^ 2 - Mod(lam1) ^ 2)
    G <- (Conj(lam0) / omega0 ^ 2) * (1 + b0) * (1 - p_star) + 
        (Conj(lam1) / omega1 ^ 2) * (1 + b1) * p_star
    return(diag(G, length(p_star)))
}

M_beta_noncir <- function(XtYre, XtYim, XtX.cv, Dv, Gv_re, Gv_im) {
    # d_star: dim = q = qq / 2
    # n <- nrow(X)
    # qq <- ncol(X)
    # return a complex-valued reg coef
    
    E <- XtX.cv + 2 * Dv
    if (dim(XtX.cv)[1] == 1) {
        b <- (2 * Gv_im * XtYre - (E - 2 * Gv_re) * XtYim) / 
             ((2 * Gv_im) ^ 2 - E ^ 2 + (2 * Gv_re) ^ 2)
        a <- (-2 * Gv_im * XtYim + (E + 2 * Gv_re) * XtYre) / 
            (E ^ 2 - (2 * Gv_re) ^ 2 - (2 * Gv_im) ^ 2)
    } else {
        denom_b <- (2 * Gv_im) %*% (2 * Gv_im) - E %*% E + (2 * Gv_re) %*% (2 * Gv_re)
        b <- chol2inv(chol(denom_b)) %*% (2 * Gv_im %*% XtYre - (E - 2 * Gv_re) %*% XtYim)
        denom_a <- E %*% E - (2 * Gv_re) %*% (2 * Gv_re) - (2 * Gv_im) %*% (2 * Gv_im)
        a <- chol2inv(chol(denom_a)) %*% ((E + 2 * Gv_re) %*% XtYre - 2 * Gv_im %*% XtYim)
    }
    return(a + b * 1i)
}



M_beta_cp_spikecir <- function(XtY, XtX, Bspike_v, Gspike_v) {
    # d_star: dim = q = qq / 2
    # n <- nrow(X)
    # qq <- ncol(X)
    if (ncol(X) > nrow(X)) {
        X_star = t(t(X) / (tausq <- rep(as.numeric(d_star), 2)))
        A = solve(diag(n) + X_star %*% t(X), tol = 1e-400)
        result = ((1 / tausq) * diag(length(tausq)) - t(X_star) %*% A %*% X_star)
        result = result %*% XtY
        return(result)
    } else {
        Psi <- chol2inv(chol(XtX + Bspike_v))
        return(Psi %*% XtY)
    }
}