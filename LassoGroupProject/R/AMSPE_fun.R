#' Function for calculating the Average Mean Squared Predicion Error (AMSPE) of an estimator
#' 
#' This function calculates the Average Mean Squared Prediction Error (AMSPE) of an estimator
#' from a Monte Carlo simulations study. 
#' @param modelFit A matrix of fitting values of dimensions number of observations x 
#' number of simulations
#' @param truth The true values of the data generating process, i.e., \eqn{x_i^T \beta}
#' @return A numeric value
#' @keywords LassoGroupProject
#' @export
AMSPE_fun <- function(modelFit, truth) {
    n <- length(truth)
    nSim <- ncol(modelFit)
    truthMat <- matrix(rep(truth, nSim), n, nSim)
    errorMat <- modelFit - truthMat
    sqErrorMat <- errorMat^2
    meanSqError <- rowMeans(sqErrorMat)
    aMeanSqError <- mean(meanSqError)
    return(aMeanSqError)
}
