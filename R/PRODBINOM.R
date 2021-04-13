#' @title Product of Independent Binomail Distributions
#' 
#' @description Density, distribution, quantile, and random generation functions for products of independent binomial distributions.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param size vector of trials (zero or more).
#' @param prob vector of probability of success on each trial.
#' @param log,log.p logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tails logical; if \code{TRUE} (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @export 
dprodbinom <- function(x, size, prob, log = FALSE) {
    ##
    if (length(prob) > 0) {
        stop("The length of 'prob' has to be larger than 1.")
    }
    
    if (length(size) > 0) {
        stop("The length of 'size' has to be larger than 1.")
    }
    
    if (length(prob) != length(size)) {
        stop("The length of 'prob' has to be equal to the length of 'size'.")
    }
    
    if (length(prob) > 2) {
        stop("Only allows for the product of two independent binomial distributions.")
    }
    
    ##
    if (any(size < 0)) {
        stop("'size' has to be larger than, or equal to, 0.")
    }
    
    if (any(prob < 0) || any(prob > 1)) {
        stop("'prob' has to be between 0 and 1.")
    }
    
    ## 
    if (any(x < 0)) {
        stop("'x' has to be larger than, or equal to, 0.")
    }
    
    ##
    p <- sapply(x, product_set_probability, N = size, p = prob)
    
    if (log) {
        p <- log(p)
    }
    
    return(p)
}

#' @rdname dprodbinom
#' 
#' @export
pprodbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
    ##
    if (length(prob) > 0) {
        stop("The length of 'prob' has to be larger than 1.")
    }
    
    if (length(size) > 0) {
        stop("The length of 'size' has to be larger than 1.")
    }
    
    if (length(prob) != length(size)) {
        stop("The length of 'prob' has to be equal to the length of 'size'.")
    }
    
    if (length(prob) > 2) {
        stop("Only allows for the product of two independent binomial distributions.")
    }
    
    ##
    if (any(size < 0)) {
        stop("'size' has to be larger than, or equal to, 0.")
    }
    
    if (any(prob < 0) || any(prob > 1)) {
        stop("'prob' has to be between 0 and 1.")
    }
    
    ## 
    if (any(q < 0)) {
        stop("'q' has to be larger than, or equal to, 0.")
    }
    
    ##
    p <- sapply(q, product_set_cumulative_probability, N = size, p = prob)
    if (!lower.tail) {
        p <- 1 - p
    }
    
    if (log.p) {
        p <- log(p)
    }
    
    return(p)
}

#' @rdname dprodbinom
#' 
#' @export
qprodbinom <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE) {
    ##
    if (length(prob) > 0) {
        stop("The length of 'prob' has to be larger than 1.")
    }
    
    if (length(size) > 0) {
        stop("The length of 'size' has to be larger than 1.")
    }
    
    if (length(prob) != length(size)) {
        stop("The length of 'prob' has to be equal to the length of 'size'.")
    }
    
    if (length(prob) > 2) {
        stop("Only allows for the product of two independent binomial distributions.")
    }
    
    ##
    if (any(size < 0)) {
        stop("'size' has to be larger than, or equal to, 0.")
    }
    
    if (any(prob < 0) || any(prob > 1)) {
        stop("'prob' has to be between 0 and 1.")
    }
    
    ##
    if (any(p < 0) || any(p > 1)) {
        stop("'p' has to be between 0 and 1.")
    }
    
    ## 
    if (log.p) {
        p <- exp(p)
    }
    
    if (!lower.tail) {
        p <- 1 - p
    }
    
    q <- sapply(p, product_set_quantile, N = size, p = prob)
    return(q)
}

#' @rdname dprodbinom
#' 
#' @export
rprodbinom <- function(n, size, prob) {
    ##
    if (length(prob) > 0) {
        stop("The length of 'prob' has to be larger than 1.")
    }
    
    if (length(size) > 0) {
        stop("The length of 'size' has to be larger than 1.")
    }
    
    if (length(prob) != length(size)) {
        stop("The length of 'prob' has to be equal to the length of 'size'.")
    }
    
    if (length(prob) > 2) {
        stop("Only allows for the product of two independent binomial distributions.")
    }
    
    ##
    if (any(size < 0)) {
        stop("'size' has to be larger than, or equal to, 0.")
    }
    
    if (any(prob < 0) || any(prob > 1)) {
        stop("'prob' has to be between 0 and 1.")
    }
    
    ##
    if (length(n) > 1) {
        n <- length(n)
    }
    
    if ((n < 1) || ((n - floor(n)) > 1e-6)) {
        stop("'n' has to be a positive integer.")
    }
    
    ## 
    K <- length(prob) 
    p <- matrix(NA, nrow = n, ncol = K)
    for (k in seq_len(K)) {
        p[, k] <- rbinom(n, size[k], prob[k])
    }
    
    p <- apply(p, 1, prod)
    return(p)
}

