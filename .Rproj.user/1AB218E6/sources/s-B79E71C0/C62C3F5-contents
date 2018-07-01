# E-step_beta_binom

library(mvtnorm)
library(LaplacesDemon)


E_beta_binom_cp_noncir <- function(beta_k, sigma_k, v0, v1, Sig_beta, theta_k) {
    qq <- length(beta_k)
    stars <- matrix(0, qq / 2, 2)
    
    # Compute p_star
    # Note: need to revise it if dealing with qq > 2
    dens1 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = sigma_k ^ 2 * v0 * Sig_beta, log = TRUE)
    dens0 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = sigma_k ^ 2 * v1 * Sig_beta, log = TRUE)
    
    stars[, 2] <- 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))
    
    # Compute d_star
    stars[, 1] <- stars[, 2] / v1 + (1 - stars[, 2]) / v0  
    return(stars)
}


E_beta_binom_cp_noncir <- function(beta_k, sigma_k, v0, v1, cor.ri, theta_k) {
    qq <- length(beta_k)
    stars <- matrix(0, qq / 2, 2)
    
    # Compute p_star
    # Note: need to revise it if dealing with qq > 2
    dens1 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = sigma_k ^ 2 * v1 * create_Dmat(cor.ri, qq / 2), 
                     log = TRUE)
    dens0 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = sigma_k ^ 2 * v0 * create_Dmat(cor.ri, qq / 2),
                     log = TRUE)
    
    stars[, 2] <- 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))
    
    # Compute d_star
    stars[, 1] <- stars[, 2] / v1 + (1 - stars[, 2]) / v0  
    return(stars)
}


library(mvnfast)
E_beta_binom_noncir <- function(beta_k, sigma_k, theta_k, omega0, omega1, 
                                   lam0, lam1) {

    # return p_stars which is the expectation of the idicator variable
    qq <- length(beta_k)
    mu <- rep(0, qq)
    
    # Compute p_star
    realcov0 <- cplx_to_real_cov(Omega = omega0, Lam = lam0)
    realcov1 <- cplx_to_real_cov(Omega = omega1, Lam = lam1)
    
    dens1 <- dmvnorm(beta_k, mean = mu, sigma = sigma_k ^ 2 * realcov1, log = TRUE)
    dens0 <- dmvnorm(beta_k, mean = mu, sigma = sigma_k ^ 2 * realcov0, log = TRUE)
    p_stars <- 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))

    return(p_stars)
}


E_beta_binom_cp_spikecir <- function(beta_k, sigma_k, v0, v1, cor.ri, theta_k) {
    qq <- length(beta_k)
    stars <- matrix(0, qq / 2, 3)
    
    # Compute p_star
    # Note: need to revise it if dealing with qq > 2
    dens1 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = sigma_k ^ 2 * v1 * create_Dmat(cor.ri, qq / 2), 
                     log = TRUE)
    dens0 <- dmvnorm(beta_k, mean = rep(0, qq), 
                     sigma = diag(sigma_k ^ 2 * v0, qq),
                     log = TRUE)
    
    stars[, 1] <- 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))
    
    # Compute diagonal d_star
    stars[, 2] <- stars[, 1] / (v1 * (1 - cor.ri ^ 2)) + (1 - stars[, 1]) / v0  
    
    # Compute off diagonal expectation
    stars[, 3] <- stars[, 1] / v1 * (-cor.ri / (1 - cor.ri ^ 2))
    
    return(stars)
}

E_beta_binom_spikecir <- function(beta_k, sigma_k, theta_k, omega0, omega1, lam0, lam1) {
    qq <- length(beta_k)
    sig2 <- sigma_k ^ 2
    mu <- rep(0, qq)
    # Compute p_star
    realcov0 <- cplx_to_real_cov(Omega = omega0, Lam = lam0)
    realcov1 <- cplx_to_real_cov(Omega = omega1, Lam = lam1)
    
        dens1 <- dmvnorm(beta_k, mean = mu, sigma = sig2 * realcov1, log = TRUE)
        dens0 <- dmvnorm(beta_k, mean = mu, sigma = sig2 * realcov0, log = TRUE)
    
    p_stars <- 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))
    
    return(p_stars)
}













