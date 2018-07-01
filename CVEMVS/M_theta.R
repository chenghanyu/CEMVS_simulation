# M-step for the probability of success of the binary variable (theta, pii)

M_theta <- function(post, a, b) {
    result = sum(post) + a - 1
    result = result / (b + a + length(post) - 2)
    return(result)
}

M_theta_cp <- function(postprob, a_th, b_th) {
    result <- (sum(postprob) + a_th - 1) / (b_th + a_th + length(postprob) - 2)
    return(result)
}