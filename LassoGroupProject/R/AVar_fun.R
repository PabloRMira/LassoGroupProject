#' Compute the Average Variance of an estimator
#' 
#' Computes the Average Variance (AV) of an estimator
#' in a Monte Carlo simulations study. 
#' @param modelFit A matrix of fitting values of dimensions (number of observations x 
#' number of simulations)
#' @return A numeric value
#' @keywords LassoGroupProject
#' @export
AVar_fun <- function(modelFit) {
    nSim <- ncol(modelFit)
    n <- nrow(modelFit)    
    expectedFit <- rowMeans(modelFit)
    expFitMat <- matrix(rep(expectedFit, nSim), n, nSim)
    sqDev <- (modelFit - expFitMat)^2
    indVar <- rowMeans(sqDev)
    aveVar <- mean(indVar)
    return(aveVar)
}
