#' Calculate the Average Squared Bias of an estimator
#' 
#' It calculates the Average Squared Bias of an estimator
#' from a Monte Carlo simulations study. 
#' @param modelFit A matrix of fitting values of dimensions number of observations x 
#' number of simulations
#' @param truth The true values of the data generating process, i.e., \eqn{x_i^T \beta}
#' @return A numeric value
#' @keywords LassoGroupProject
#' @export
ABias_fun <- function(modelFit, truth) {
    expectedFit <- rowMeans(modelFit)
    sqBias <- (expectedFit - truth)^2
    aveSqBias <- mean(sqBias)
    return(aveSqBias)
}
