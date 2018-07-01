# log prior of the binary variable

log_prior <- function (Gamma, a, b, NN) {
    x <- sum(Gamma) + a
    y <- NN - sum(Gamma) + b
    prior_val <- base::lbeta(x, y) - (base::lbeta(a, b))
    
    # Stirling approximation
    if(prior_val == "Inf" || prior_val == "-Inf"){
        prior_val = 0.5 * log(2 * pi) + (x - 0.5) * log(x) + 
            (y - 0.5) * log(y) - (x + y - 0.5) * log(x + y)
    }
    return(prior_val)
}