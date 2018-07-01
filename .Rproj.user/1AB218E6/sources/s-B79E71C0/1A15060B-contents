
# Circular EMVS one theta one sig main function
EMmod_theta_sig <- function(X, Y, v0s, v1, epsilon, theta, a, b, thre) {
    # X: n x q matrix 
    # Y: n x N x N Array
    # v0: vector of length L
    # v1: scalar
    # beta_init: q x 1 vector
    # sigma_init: scalar
    # epsilon: scalar
    # theta: scalar
    # a: scalar
    # b: scalar
    
    if (missing(theta)){theta = 0.5}
    
    if (class(X) != "matrix") {
        tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("X must be a matrix or able to be coerced to a matrix")
    }
    
    if (class(Y) != "array") {
        tmp <- try(Y <- as.array(Y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("Y must be an array or able to be coerced to an array")
    }
    
    if (sum(v1<v0s)>0)
        stop("v1 has to be larger than v0")
    
    if (length(a) > 1) 
        stop("a has to be a scalar")
    
    if (length(b) > 1) 
        stop("b has to be a scalar")
    
    if (length(v1) > 1) 
        stop("v1 has to be a scalar")
    
    if (length(theta) > 1) 
        stop("theta has to be a scalar")
    
    q <- ncol(X)
    n <- nrow(X)
    
    if(n!=dim(Y)[1])
        stop("length of Y for each voxel has to be nrow(X)")
    
    # n <- dim(Y)[1] # 2 * length of time series 
    N <- dim(Y)[2] # matrix size
    # q <- ncol(X) # number of reg coeff,  number of cplx = 2 * number of mag 
    XtX <- t(X) %*% X
    Xt <- t(X)
    L <- length(v0s)
    # Arrays for results
    sigmas = rep(0, L)
    betas = array(0, c(L, q, N, N))
    beta_k = array(0, c(q, N, N))
    beta_new = array(0, c(q, N, N))
    posts = array(0, c(L, q, N, N)) 
    p_star = array(0, c(q, N, N))
    d_star = array(0, c(q, N, N))
    log_gfcn = rep(0, L)
    counts = rep(0, L)
    thetas = rep(0, L)
    
    Beta_init = array(0, c(q, N, N))
    # Sigma_init = 1
    
    for (m in 1:L) {
        v0 = v0s[m]
        print(paste("v0 =", v0))
        eps = epsilon + 1
        count <- 1
        # p_star_sum = 0
        # Initialization of parameters
        theta_k = 0.5
        sigma_k = 1
        for (j in 1:N) {
            for (i in 1:N) {
                # initial values: Start with LSE
                lmout <- lm(Y[, i, j] ~ 0 + X)
                Beta_init[, i, j] <- lmout$coeff
                beta_k[, i, j] = Beta_init[, i, j]
                beta_new[, i, j] = beta_k[, i, j]
            }
        }
        sigma_new = sigma_k
        theta_new = theta_k
        while(eps > epsilon) {
            if (count == 2000) break
            
            for (j in 1:N) {
                for (i in 1:N) {
                    # ******** E-step ************ #
                    # Update inclusion probabilities
                    # E_step contains p_star and d_star
                    E_step = E_beta_binom(beta_k[, i, j], sigma_k, 
                                          v0, v1, theta_k)
                    d_star[, i, j] = E_step[, 1]
                    p_star[, i, j] = E_step[, 2]
                    
                    # ******** M-step for beta ************ #
                    XtY = Xt %*% Y[, i, j]
                    beta_k[, i, j] = M_beta(XtY, X, Y[, i, j], XtX, 
                                            d_star[, i, j])
                }
            }
            # ******** M-step for theta and sig ************ #
            theta_k <- M_theta(p_star, a, b)
            sigma_k <- M_sigma_sig(Y, X, beta_k, d_star, 1, 1) #check
            eps = sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            beta_new = beta_k
            sigma_new = sigma_k
            theta_new = theta_k
            count = count + 1
            if (count %% 2 == 0) print(paste("Epsilon =", round(eps, 4), 
                                             " count", count))
            if (eps < epsilon) print(paste("final count =", count))
        }
        
        # Store values:
        posts[m, , , ] = p_star
        betas[m, , , ] = beta_new
        sigmas[m]  = sigma_new
        thetas[m] = theta_new
        counts[m] = count
        Gamma = p_star > thre
        log_gfcn[m] = log_g_sig(Gamma, X, Y, eta = 1, lambda = 1, v0 = 0, 
                                v1, a, b)
    }
    return(list(pprob = posts, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts))
}

EM_theta_sig <- function(X, Y, v0s, v1, epsilon, theta, a, b, thre) {
    # X: n x qq matrix 
    # Y: n x N x N Array
    # v0: vector of length L
    # v1: scalar
    # beta_init: qq x 1 vector
    # sigma_init: scalar
    # epsilon: scalar
    # theta: scalar
    # a: scalar
    # b: scalar
    # va_g: scalar
    
    if (missing(theta)){theta = 0.5}
    
    if (class(X) != "matrix") {
        tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("X must be a matrix or able to be coerced to a matrix")
    }
    
    if (class(Y) != "array") {
        tmp <- try(Y <- as.array(Y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("Y must be an array or able to be coerced to an array")
    }
    
    if (sum(v1<v0s) > 0)
        stop("v1 has to be larger than v0")
    
    if (length(a) > 1) 
        stop("a has to be a scalar")
    
    if (length(b) > 1) 
        stop("b has to be a scalar")
    
    if (length(v1) > 1) 
        stop("v1 has to be a scalar")
    
    if (length(theta) > 1) 
        stop("theta has to be a scalar")
    
    qq <- ncol(X)
    n <- nrow(X)
    q <- qq / 2
    
    if(n != dim(Y)[1])
        stop("length of Y for each voxel has to be nrow(X)")
    
    N <- dim(Y)[2] # matrix size
    XtX <- t(X) %*% X
    Xt <- t(X)
    L <- length(v0s)
    # Arrays for results
    sigmas = rep(0, L)
    betas = array(0, c(L, qq, N, N))
    beta_k = array(0, c(qq, N, N))
    beta_new = array(0, c(qq, N, N))
    posts = array(0, c(L, q, N, N)) 
    p_star = array(0, c(q, N, N))
    d_star = array(0, c(q, N, N))
    log_gfcn = rep(0, L)
    counts = rep(0, L)
    thetas = rep(0, L)
    
    Beta_init = array(0, c(qq, N, N))
    Sigma_init = array(0, c(N, N))
    
    # v1 <- rep(v1, qq)
    
    for (m in 1:L) {
        v0 = v0s[m]
        print(paste("v0 =", v0))
        eps = epsilon + 1
        count <- 1
        # p_star_sum = 0
        # Initialization of parameters
        theta_k = 0.5
        sigma_k = 1
        for (j in 1:N) {
            for (i in 1:N) {
                # initial values: Start with LSE
                lmout <- lm(Y[, i, j] ~ 0 + X)
                Beta_init[, i, j] <- lmout$coeff
                beta_k[, i, j] = Beta_init[, i, j]
                beta_new[, i, j] = beta_k[, i, j]
            }
        }
        sigma_new = sigma_k
        theta_new = theta_k
        while(eps > epsilon) {
            if (count == 2000) break
            
            for (j in 1:N) {
                for (i in 1:N) {
                    # ******** E-step ************ #
                    # Update inclusion probabilities
                    # E_step contains p_star and d_star
                    E_step = E_beta_binom_cp(beta_k[, i, j], sigma_k, 
                                             v0, v1, theta_k)
                    d_star[, i, j] = E_step[, 1]
                    p_star[, i, j] = E_step[, 2]
                    # p_star_sum = p_star_sum + sum(p_star[, i, j])
                    
                    # ******** M-step for beta ************ #
                    XtY = Xt %*% Y[, i, j]
                    beta_k[, i, j] = M_beta_cp(XtY, X, Y[, i, j], XtX, 
                                               d_star[, i, j])
                }
            }
            # ******** M-step for theta and sigma ************ #
            theta_k <- M_theta_cp(p_star, a, b)
            sigma_k <- M_sigma_sig_cp(Y, X, beta_k, d_star, 1, 1) #check
            eps = sum((beta_new - beta_k) ^ 2) + sum((sigma_new - sigma_k) ^ 2) + 
                sum((theta_new - theta_k) ^ 2)
            beta_new = beta_k
            sigma_new = sigma_k
            theta_new = theta_k
            count = count + 1
            if (count %% 2 == 0) print(paste("Epsilon =", round(eps, 4), 
                                             " count", count))
            if (eps < epsilon) print(paste("final count =", count))
        }
        
        # Store values:
        posts[m, , , ] = p_star
        betas[m, , , ] = beta_new
        sigmas[m] = sigma_new
        thetas[m] = theta_new
        counts[m] = count
        Gamma = (p_star > thre)
        log_gfcn[m] = log_g_sig_cp(Gamma, X, Y, eta = 1, lambda = 1, v0 = 0, 
                                   v1, a, b)
    }
    return(list(pprob = posts, betas = betas, sigmas = sigmas, 
                thetas = thetas, log_g_function = log_gfcn, v0 = v0s,
                v1 = v1, counts = counts))
}

# Other functions used in Analysis
summaryplot_sig <- function(sig, result, type, v0s) {
    for (L in 1:length(v0s)) {
        image.plot(result$pprob[L, 1, , ], axes=FALSE,
                   zlim = c(0, 1), main = paste(type, "v0=", round(v0s[L], 4),
                                                "\n L=", L))
        print(paste("v0=", v0s[L]))
        print(paste("L=", L, "sig =", result$sigmas[L], "theta=", result$thetas[L]))
    }
}

# Function to plot an activation map
imageplot <- function(aa, posprob, main) {
    for (i in aa) {
        image(posprob > i, main = paste(main), axes = FALSE,
              col = c("navy", "red3") )
    }
}


# Function to plot activation maps and strength maps for CP and MO models
act_str_plot <- function (datatype, snr, aa, islabel) {
    # datatype: a number among 1 to 12
    par(mfrow = c(1,1))
    par(mar = c(.5, .5, 2, .5))
    L = which(v0s_mod == EMVSbest(em_theta_sig[[datatype]])$v0)
    imageplot(aa, em_theta_sig[[datatype]]$pprob[L, 1, , ], 
              main = paste("CP Activation SNR =", snr))
    Lmod  = which(v0s_mod == EMVSbest(emmod_theta_sig[[datatype]])$v0)
    imageplot(aa, emmod_theta_sig[[datatype]]$pprob[Lmod, 1, , ], 
              main = paste("MO Activation SNR =", snr))
    
    par(mar = c(.5, .5, 2, .5))
    betas = em_theta_sig[[datatype]]$betas[L, , , ]
    strength_beta = sqrt(betas[1, , ] ^ 2 + betas[2, , ] ^ 2)
    if (islabel) {
        # check strength of voxels classified as active
        strength_beta[em_theta_sig[[datatype]]$pprob[L, 1, , ] < aa] <- 0
    } else {
        # check strength of true active voxels
        strength_beta[final == 0] <- 0 
    }
    image.plot(strength_beta, axes=FALSE, cex.main = 2, legend.cex = 1.5,
               zlim = c(0, 1), main = paste("CP Strength SNR =", snr))
    betasmod = emmod_theta_sig[[datatype]]$betas[Lmod, , , ]
    if (islabel) {
        # check strength of voxels classified as active
        betasmod[emmod_theta_sig[[datatype]]$pprob[Lmod, 1, , ] < aa] <- 0
    } else {
        # check strength of true active voxels
        betasmod[final == 0] <- 0
    }
    image.plot(abs(betasmod), axes=FALSE, cex.main = 2, legend.cex = 1.5,
               zlim = c(0, 1), main = paste("MO Strength SNR =", snr))
}


# Function to compute true positives and false positives
tpfpDatafcn <- function(cutoff, prob, truedata) {
    tpfp <- function(cutoff, prob, truedata) {
        # compute true positive rate and false positive rate
        for (i in 1:length(cutoff)) {
            tbl <- table(prob > cutoff[i], truedata > 0)
            marg <- margin.table(tbl, 2)
            if ("FALSE" %in% rownames(tbl) && "TRUE" %in% rownames(tbl)) {
                # true positive
                tp <- tbl[2, 2] / marg[2]
                # false positive
                fp <- tbl[2, 1] / marg[1]
            } else if (!("FALSE" %in% rownames(tbl))) {
                # true positive
                tp <- tbl[1, 2] / marg[2]
                # false positive
                fp <- tbl[1, 1] / marg[1]
            } else if  (!("TRUE" %in% rownames(tbl))) {
                # true positive
                tp <- 0
                # false positive
                fp <- 0
            }
            rate <- c(tp, fp)
        }
        tbl <- table(prob > cutoff, truedata)
        marg <- margin.table(tbl, 2)
        
        names(rate) <- c('TruePos','FalsePos')
        return(rate)
    }
    tpfpData <- as.data.frame(
        do.call(rbind, lapply(cutoff, tpfp, prob, truedata)))
    tpfpData$threshold <- cutoff
    return(tpfpData)
}

# Function used for lasso
tpfpLasso <- function(betatf, truedata) {
    # truedata: a true/false matrix
    # compute true positive rate and false positive rate
    tbl <- table(betatf, truedata)
    marg.col <- margin.table(tbl, 2)
    marg.row <- margin.table(tbl, 1)
    if ("FALSE" %in% rownames(tbl) && "TRUE" %in% rownames(tbl)) {
        # true positive (recall)
        tp <- tbl[2, 2] / marg.col[2]
        # false positive
        fp <- tbl[2, 1] / marg.col[1]
        # precision
        precision <- tbl[2, 2] / marg.row[2]
        # accuracy
        acc = (tbl[1, 1] + tbl[2, 2]) / sum(tbl)
    } else if (!("FALSE" %in% rownames(tbl))) {
        # true positive
        tp <- tbl[1, 2] / marg.col[2]
        # false positive
        fp <- tbl[1, 1] / marg.col[1]
        # precision
        precision <- tbl[1, 2] / marg.row[1]
        # accuracy
        acc = tbl[1, 2] / sum(tbl)
    } else if  (!("TRUE" %in% rownames(tbl))) {
        # true positive
        tp <- 0
        # false positive
        fp <- 0
        # precision
        precision <- 0
        # accuracy
        acc = tbl[1, 1] / sum(tbl)
    }
    rate <- c(tp, 1 - tp, fp, 1 - fp, precision, 1 - precision, acc)
    names(rate) <- c('TPR/sensitivity/recall', "FNR=1-TPR",
                     'FPR', 
                     'TNR=1-FPR(specificity)',
                     "Precision(PPV)", "FDR=1-PPV", "Accuracy(ACC)")
    return(rate)
}

# Function to plot ROC curve (from EMbeta.R script)
ROC <- function(tpfpDataList, sig2, r, c) {
    plot(tpfpDataList[[1]]$FalsePos, tpfpDataList[[1]]$TruePos, 
         type = 's', ylab = "True Positive Rate", xlab = "False Positive Rate", 
         col = 1, lwd = 2, 
         main = paste("ROC, sig2 =", sig2, ", ", "(r, c) = (", r, ",", c, ")", 
                      sep = ""),
         xlim = c(0, 1), ylim = c(0, 1))
    lines(tpfpDataList[[2]]$FalsePos, tpfpDataList[[2]]$TruePos, type = 's', 
          col = 2, lty = 2, lwd = 2)
    lines(tpfpDataList[[3]]$FalsePos, tpfpDataList[[3]]$TruePos, type = 's', 
          col = 3, lty = 3, lwd = 2)
    lines(tpfpDataList[[4]]$FalsePos, tpfpDataList[[4]]$TruePos, type = 's', 
          col = 4, lty = 4, lwd = 2)
    legend('bottomright', c('Cplx', 'Re', 'Im', 'Mod'), 
           col = 1:4, lwd = rep(2, 4), lty = 1:4, cex = .9, bty = 'n')
}

# functions for plotting brain image activations 
img_fmri <- function(data.mean.list, output, thre, v0s, v1, iscplx, isstr, title) {
    # par(mfrow = c(1, 1))
    opt_v0 = EMVSbest(output)$v0
    opt_v0_idx = which(v0s == EMVSbest(output)$v0)
    if (missing(title)) {
        title = paste(ifelse(iscplx, "EM CV v0 =", "EM MO v0 ="), round(opt_v0, 4))
    }
    image(data.mean.list$Mod, col = gray((0:96) / 96), axes = FALSE, main = title)
    posprob <- output$pprob[opt_v0_idx, 1, , ]
    for (i in 1:96) {
        ll <- length(which(posprob[i, ] > thre))
        idx <- which(posprob[i, ] > thre)
        if (ll > 0) {
            if(isstr) {
                par(mar = c(.5, .5, 2, 5))
                danga <- output$betas[opt_v0_idx, , , ]
                # strength = sqrt(re^2 + im^2)
                if (iscplx) {
                    danga <- sqrt(danga[1, , ] ^ 2 + danga[2, , ] ^ 2) 
                    rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
                } else {
                    rbPal <- colorRampPalette(c('darkblue', 'purple', 'red2','orange2'))
                }
                gaact = danga[which(posprob > thre)]
                brk = 10
                #                 Col1 <- rbPal(brk)[as.numeric(cut(seq(-0.5, 1, length = 9216), 
                #                                                   breaks = brk))]
                Col1 <- rbPal(brk)[as.numeric(cut(as.vector(gaact), breaks = brk))]
                Color = matrix(0, nc = 96, nr = 96)
                Color[which(posprob > thre)] = Col1
                points(rep((i - .5)/96, ll), (which(posprob[i, ] > thre) - .5)/96, 
                       col = Color[i, idx], pch = 15, cex = .8)
                mini = min(gaact * 0.95)
                maxi = max(gaact * 1.05)
                subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
                        x = 1.01, y = .5, size = c(1, 6.7))
            } else {
                # par(omi = c(0, 0, 0, 0))
                # par(mar = c(.5,.5, 2, .5))
                points(rep((i - .5)/96, ll), (which(posprob[i, ] > thre) - .5)/96, 
                       col = "red", pch = 15, cex = 0.8)
            }
        }
    }
}

# functions for plotting the optimal brain image activations 
img_fmri <- function(braindata = data.mean.list$Mod,
                     output, thre, v0s, v1, iscplx, isstr, title,
                     isrotate = FALSE, isbar = FALSE, iscut = FALSE,
                     cex.pt = 0.69, bar.loc.x = 0.99, bar.loc.y = 0.5,
                     bar.size = c(0.5, 6), cex.main = 1.2,
                     cut.row.idx = c(1:20, 77:96),
                     cut.col.idx = c(1:15, 82:96)){
    scale_fcn <- function(x, a, b) {
        (b - a) * (x - min(x)) / (max(x) - min(x))
    }
    # par(mfrow = c(1, 1))
    if (isbar) {
        # par(mar = c(.85, .5, 2, 5))
        par(mar = c(.1, .1, 1.5, 5))
    } else {
        par(mar = c(.1, .1, 1.5, .1))
        # par(mar = c(.5, .5, 2, .5))
    }
    opt_v0 <- EMVSbest(output)$v0
    opt_v0_idx <- which(v0s == EMVSbest(output)$v0)
    if (missing(title)) {
        title <- paste(ifelse(iscplx, "EM CV v0 =", "EM MO v0 ="), round(opt_v0, 4))
    }
    if (isrotate) {
        posprob <- rotate(output$pprob[opt_v0_idx, 1, , ])
    } else {
        posprob <- output$pprob[opt_v0_idx, 1, , ]
    }
    danga <- output$betas[opt_v0_idx, , , ]
    brk <- 120
    for(k in 1:length(thre)) {
        #         image(braindata, col = gray((0:96) / 96), axes = FALSE, 
        #               main = paste(title, "thre =", thre[k]))
        #                     # strength = sqrt(re^2 + im^2)
        if (iscplx) {
            danga <- sqrt(danga[1, , ] ^ 2 + danga[2, , ] ^ 2) 
            # rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
            rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        } else {
            rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        }
        if(iscut) {
            # cut.row.idx <- c(1:20, 77:96)
            # cut.col.idx <- c(1:15, 82:96)
            posprob <- posprob[-cut.row.idx, -cut.col.idx]
            # braindata <- braindata[-cut.row.idx, -cut.col.idx]
            danga <- danga[-cut.row.idx, -cut.col.idx]
        }
        image(braindata, col = gray((0:96) / 96), axes = FALSE, 
              main = paste(title), cex.main = cex.main)
        for (i in 1:nrow(posprob)) {
            ll <- length(which(posprob[i, ] > thre[k]))
            idx <- which(posprob[i, ] > thre[k])
            if (ll > 0) {
                if(isstr) {
                    danga <- scale_fcn(danga, 0, 0.5)
                    gaact <- danga[which(posprob > thre)]
                    #                 Col1 <- rbPal(brk)[as.numeric(cut(seq(-0.5, 1, length = 9216), 
                    #                                                   breaks = brk))]
                    Col1 <- rbPal(brk)[as.numeric(cut(as.vector(gaact), breaks = brk))]
                    Color <- matrix(0, nrow = nrow(braindata), ncol = ncol(braindata))
                    Color[which(posprob > thre[k])] <- Col1
                    # Color[which(true.strength.scale5 > 0)] <- Col1
                    points(rep(i/nrow(braindata), ll), 
                           (which(posprob[i, ] > thre[k]))/ncol(braindata), 
                           col = Color[i, idx], pch = 15, cex = cex.pt)
                    #                     points(rep((i - .5)/96, ll), (which(true.strength.scale5[i, ] > 0) - .5)/96, 
                    #                            col = Color[i, idx], pch = 15, cex = .69)
                    mini <- min(gaact * 0.95)
                    maxi <- max(gaact * 1.05)
                } else {
                    # par(omi = c(0, 0, 0, 0))
                    # par(mar = c(.5,.5, 2, .5))
                    points(rep(i / nrow(braindata), ll), 
                           (which(posprob[i, ] > thre[k]))/ncol(braindata), 
                           col = "red", pch = 15, cex = cex.pt)
                }
            }
        }
        if (isstr) {
            #             subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
            #                     x = .99, y = .5, size = c(1, 6))
            if(isbar) {
                subplot(color.bar(rbPal(brk), min = 0, max = 0.5, cex.axis = 1), 
                        x = bar.loc.x, y = bar.loc.y, 
                        size = bar.size)   
            }
        }
    }
}

# functions for plotting the optimal brain image activations 
img_fmri <- function(braindata = data.mean.list$Mod, output,
                     thre, v0s, v1, iscplx, isstr, title,
                     isrotate = FALSE, isbar = FALSE, iscut = FALSE,
                     cex.pt = 0.69, bar.loc.x = 0.99, bar.loc.y = 0.5,
                     bar.size = c(0.5, 6), cex.main = 1.2,
                     cut.row.idx = c(1:20, 77:96),
                     cut.col.idx = c(1:15, 82:96),
                     cpp = FALSE){
    scale_fcn <- function(x, a, b) {
        (b - a) * (x - min(x)) / (max(x) - min(x))
    }
    # par(mfrow = c(1, 1))
    if (isbar) {
        # par(mar = c(.85, .5, 2, 5))
        par(mar = c(.1, .1, 1.5, 5))
    } else {
        par(mar = c(.1, .1, 1.5, .1))
        # par(mar = c(.5, .5, 2, .5))
    }
    opt_v0 <- EMVSbest(output)$v0
    opt_v0_idx <- which(v0s == EMVSbest(output)$v0)
    if (missing(title)) {
        title <- paste(ifelse(iscplx, "EM CV v0 =", "EM MO v0 ="), round(opt_v0, 4))
    }
    
    if(cpp) {
        pprob <- matrix(output$pprob[1, , opt_v0_idx], dim(data.mean.list$Mod))
        danga <- array(output$betas[, , opt_v0_idx], dim = c(2, dim(data.mean.list$Mod)))
    } else {
        pprob <- output$pprob[opt_v0_idx, 1, , ]
        danga <- output$betas[opt_v0_idx, , , ]
    }
    if (isrotate) {
        posprob <- rotate(pprob)
    } else {
        posprob <- pprob
    }

    brk <- 120
    for(k in 1:length(thre)) {
        #         image(braindata, col = gray((0:96) / 96), axes = FALSE, 
        #               main = paste(title, "thre =", thre[k]))
        #                     # strength = sqrt(re^2 + im^2)
        if (iscplx) {
            danga <- sqrt(danga[1, , ] ^ 2 + danga[2, , ] ^ 2) 
            # rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
            rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        } else {
            rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        }
        if(iscut) {
            # cut.row.idx <- c(1:20, 77:96)
            # cut.col.idx <- c(1:15, 82:96)
            posprob <- posprob[-cut.row.idx, -cut.col.idx]
            # braindata <- braindata[-cut.row.idx, -cut.col.idx]
            danga <- danga[-cut.row.idx, -cut.col.idx]
        }
        image(braindata, col = gray((0:96) / 96), axes = FALSE, 
              main = paste(title), cex.main = cex.main)
        for (i in 1:nrow(posprob)) {
            ll <- length(which(posprob[i, ] > thre[k]))
            idx <- which(posprob[i, ] > thre[k])
            if (ll > 0) {
                if(isstr) {
                    danga <- scale_fcn(danga, 0, 0.5)
                    gaact <- danga[which(posprob > thre)]
                    #                 Col1 <- rbPal(brk)[as.numeric(cut(seq(-0.5, 1, length = 9216), 
                    #                                                   breaks = brk))]
                    Col1 <- rbPal(brk)[as.numeric(cut(as.vector(gaact), breaks = brk))]
                    Color <- matrix(0, nrow = nrow(braindata), ncol = ncol(braindata))
                    Color[which(posprob > thre[k])] <- Col1
                    # Color[which(true.strength.scale5 > 0)] <- Col1
                    points(rep((i - 0.5)/nrow(braindata), ll), 
                           (which(posprob[i, ] > thre[k])- 0.5)/ncol(braindata), 
                           col = Color[i, idx], pch = 15, cex = cex.pt)
                    #                     points(rep((i - .5)/96, ll), (which(true.strength.scale5[i, ] > 0) - .5)/96, 
                    #                            col = Color[i, idx], pch = 15, cex = .69)
                    mini <- min(gaact * 0.95)
                    maxi <- max(gaact * 1.05)
                } else {
                    # par(omi = c(0, 0, 0, 0))
                    # par(mar = c(.5,.5, 2, .5))
                    points(rep((i - 0.5)/ nrow(braindata), ll), 
                           (which(posprob[i, ] > thre[k])- 0.5)/ncol(braindata), 
                           col = "red", pch = 15, cex = cex.pt)
                }
            }
        }
        if (isstr) {
            #             subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
            #                     x = .99, y = .5, size = c(1, 6))
            if(isbar) {
                subplot(color.bar(rbPal(brk), min = 0, max = 0.5, cex.axis = 1), 
                        x = bar.loc.x, y = bar.loc.y, 
                        size = bar.size)   
            }
        }
    }
}

color.bar <- function(lut, min, max = -min, nticks = 11, 
                      ticks = seq(min, max, len = nticks), title = '',
                      cex.axis = .8) {
    scale = (length(lut) - 1) / (max - min)
    
    # dev.new(width=1.75, height=5)
    #     par(mfrow = c(1, 1))
    # par(mar = c(10,2,3,1))
    plot(c(0, 1), c(min, max), type = 'n', bty = 'n', xaxt = 'n', 
         xlab = '', yaxt = 'n', ylab = '', main = title)
    axis(4, round(ticks, 3), las = 1, cex.axis = cex.axis)
    for (i in 1:(length(lut)-1)) {
        y = (i - 1) / scale + min
        rect(0, y, 1, y + 1/scale, col = lut[i], border = NA)
    }
}
multi_brain <- function(output, thre, v0s, iscplx, title) {
    par(mfrow = c(4, 4))
    par(mar = c(.5, .5, 2, .5))
    for (j in 1:length(v0s)) {
        posprob <- output$pprob[j, 1, , ]
        danga <- output$betas[j, , , ]
        # posprob = f(posprob)
        if(missing(title)) {
            title = paste(ifelse(iscplx, "EM CV v0 =", "EM MO v0 ="), 
                           round(v0s[j], 4), title)
        }
        image(data_mean_list$Mod, col = gray((0:96) / 96), axes = FALSE, 
              main = paste(title, round(v0s[j], 4)))
        if (iscplx) {
            danga <- sqrt(danga[1, , ] ^ 2 + danga[2, , ] ^ 2) 
            rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
        } else {
            rbPal <- colorRampPalette(c('darkblue', 'purple', 
                                        'red2','orange2'))
        }
        brk = 10
        for (i in 1:96) {
            ll <- length(which(posprob[i, ] > thre))
            idx <- which(posprob[i, ] > thre)
            #         if (iscplx) {
            #             ga <- sqrt(ga[1, , ] ^ 2 + ga[2, , ] ^ 2)
            #         }
            # Color = matrix(0, nc = 96, nr = 96)
            for (s in 1:length(thre)) {
                #             gaact = ga[which(posp1 > s)] # check
                #             rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
                #             brk = 10
                #             Col1 <- rbPal(brk)[as.numeric(cut(as.vector((gaact)), 
                #                                               breaks = brk))]
                #             Color[which(posp1 > s)] = Col1
                # strength = sqrt(re^2 + im^2)
                gaact = danga[which(posprob > thre[s])]
                # rbPal <- colorRampPalette(c('red2','orange2','gold', 'yellow'))
                Col1 <- rbPal(brk)[as.numeric(cut(as.vector(gaact), breaks = brk))]
                Color = matrix(0, nc = 96, nr = 96)
                Color[which(posprob > thre[s])] = Col1
                points(rep((i - .5)/96, ll), (which(posprob[i, ] > thre[s]) - .5)/96, 
                       col = Color[i, idx], pch = 15, cex = .5)
            }
        }
    }
}

scale_fcn <- function(x, a, b) {
    (b - a) * (x - min(x)) / (max(x) - min(x))
}

# rotate() is a rotation function that gives us right plotting direction 
rotate <- function(m) t(m)[, nrow(m):1]

# Add true voxel with strength in the image
addTrueVoxelStr <- function(map, brk = 120, cex.pt = 1) {
    # if(iscut) {
    #     cut.row.idx <- c(1:20, 77:96)
    #     cut.col.idx <- c(1:15, 82:96)
    #     posprob <- posprob[-cut.row.idx, -cut.col.idx]
    #     braindata <- braindata[-cut.row.idx, -cut.col.idx]
    #     danga <- danga[-cut.row.idx, -cut.col.idx]
    # }
    # image(braindata, col = gray((0:96) / 96), axes = FALSE, 
    #       main = paste(title), cex.main = cex.main)
    Nrow <- nrow(map) 
    Ncol <- ncol(map)
    for (i in 1:Nrow) {
        ll <- length(which(map[i, ] > 0))
        idx <- which(map[i, ] > 0)
        if (ll > 0) {
            rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", 
                                        "#7F0000"))
            gaact <- map[map > 0]
            Col1 <- rbPal(brk)[as.numeric(cut(as.vector(gaact), breaks = brk))]
            Color <- matrix(0, nc = Ncol, nr = Nrow)
            Color[which(map > 0)] <- Col1
            points(rep((i - 0.5)/ Nrow, ll), 
                   (which(map[i, ] > 0) - 0.5) / Ncol, 
                   col = Color[i, idx], pch = 15, cex = cex.pt)
        }
    }
}





