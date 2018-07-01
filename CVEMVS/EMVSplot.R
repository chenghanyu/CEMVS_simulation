EMVSplot <- function(result, plot_type, omit.zeroes) {
    
    betas <- result$betas
    # intersects <- result$intersects
    log_post <- result$log_g_function
    v0s <- result$v0
    
    if (plot_type == "both"){
        par(mfrow = c(1,2))
    }
    
    if (plot_type == "both" | plot_type == "reg"){
        
        select <- apply(result$prob_inclusion, 2, 
                        function(x){as.numeric(x > 0.5)})
        
        if (omit.zeroes){betas[select == 0] = 0}
        
        matplot(v0s, betas, xlab = expression(v[0]), 
                ylab = expression(hat(beta)),
                lwd = 1, col = "grey", lty = 2, type = "l")
        matpoints(v0s, betas*select, xlab = expression(v[0]),
                  ylab = expression(hat(beta)), lwd = 1, col = 4, lty = 2, 
                  pch = 19)
        matpoints(v0s, betas*(1 - select), xlab = expression(v[0]),
                  ylab = expression(hat(beta)), lwd = 1, col = 2, lty = 2,
                  pch = 19)
        title("EMVS Regularization Plot")
        
        par(xpd = T)
        labels = paste("X", 1:ncol(betas), sep="")
        labels[select[length(v0s), ] == 0] <- ""
        text(max(v0s) * (1.1), betas[length(v0s), ], labels = labels, col = 4)
    }
    
    if (plot_type == "both" | plot_type == "gf"){
        
        plot(as.numeric(log_post) ~ v0s, pch = 4, type = "b", col = 2, lwd = 1,
             xlab = expression(v[0]), ylab = expression(log(g(gamma))))
        
        title(expression(Log(g(gamma))))}
    plot(as.numeric(log_post) ~ v0s, pch = 4, type = "b", col = 2, lwd = 1,
         xlab = expression(v[0]), ylab = expression(log(g(gamma))))
    
    title(expression(Log(g(gamma))))
}

aa <- c(0.5) 
imageplot <- function(aa, posprob, main) {
    for (i in aa) {
        image(posprob > i, main = paste(main, " thre = ", i), axes=FALSE,
              col = tim.colors()) 
    }
}

