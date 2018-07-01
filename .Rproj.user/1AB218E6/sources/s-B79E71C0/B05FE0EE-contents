# E-step for the indicator variables (circular normal case)

library(mvtnorm)
E_beta_binom <- function(beta_k, sigma_k, v0, v1, theta_k) {
    q = length(beta_k)
    stars <- matrix(0, q, 2) 
    
    # Compute p_star
    dens1 = dmvnorm(beta_k, mean = rep(0, q), 
                    sigma = diag(sigma_k ^ 2 * v1, q), log = TRUE)
    dens0 = dmvnorm(beta_k, mean = rep(0, q), 
                    sigma = diag(sigma_k ^ 2 * v0, q), log = TRUE)
    
    stars[, 2] = 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))

    # Compute d_star
    stars[, 1] =  stars[, 2] / v1 + (1 - stars[, 2]) / v0  
    return(stars)
}

E_beta_binom_cp <- function(beta_k, sigma_k, v0, v1, theta_k) {
    qq = length(beta_k)
    stars <- matrix(0, qq / 2, 2)
    
    # Compute p_star
    # Note: need to revise it if dealing with qq > 2
    dens1 = dmvnorm(beta_k, mean = rep(0, qq), 
                    sigma = diag(sigma_k ^ 2 * v1, qq), log = TRUE)
    dens0 = dmvnorm(beta_k, mean = rep(0, qq), 
                    sigma = diag(sigma_k ^ 2 * v0, qq), log = TRUE)
    
    stars[, 2] = 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1))
    
    # Compute d_star
    stars[, 1] =  stars[, 2] / v1 + (1 - stars[, 2]) / v0  
    return(stars)
}
