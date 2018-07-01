# Find the best model based on log_g function


EMVSbest <- function(result) {
    posts = result$log_g_function
    posts[is.finite(posts) == F] = NaN
    which <- which(posts == max(posts))
    logpost <- posts[which[1]]
    print("Best Model Found")
    list(log_g_function = logpost, v0 = result$v0[which[length(which)]])
}