#' Adaptive Lasso with automatic \eqn{\gamma} selection by cross-validation
#' 
#' Computes the best (\eqn{\lambda}, \eqn{\gamma}) pair for the adaptive Lasso by 
#' cross-validation.
#' @param X A matrix of regressors.
#' @param y The response variable.
#' @param gammaSeq A sequence of gamma values to cross-validate.
#' @param standardize Should the regressors be standardized?
#' @param lambdaOption Should the lambda be the minimizer of the cross-validated
#' error curve (min) or should the lambda be chosen such that its corresponding
#' cross-validated error is only 1-standard-deviation above the minimum error (1se)?
#' @param initEst Initial estimator for the adaptive Lasso. 
#' Either OLS (OLS) or ridge regression (ridge).
#' @return The tuple (\eqn{\lambda}, \eqn{\gamma}) of minimizers of the cross-validated 
#' error curve for the adaptive Lasso together with the estimated model.
#' @keywords LassoGroupProject
#' @export
AdalassoCV <- function(X, y, gammaSeq = c(0.5, 1, 2, 5),
                       standardize = FALSE,
                       lambdaOption = "min",
                       initEst = "OLS") {
  
  if (initEst == "OLS") {
    # OLS as first stage estimator
    df <- data.frame(y=y, X=X)
    olsModel <- lm(y ~ ., data = df)
    initCoefs <- as.vector(coef(olsModel))[-1]
  } else if (initEst == "ridge") {
    # Ridge regression as initial estimator
    ridgeReg <- glmnet::cv.glmnet(X, y, family="gaussian",
                                  alpha=0,
                                  standardize=standardize)
    lambda <- ridgeReg$lambda.min
    
    # Extract ridge coefficients
    initCoefs <- as.vector(coef(ridgeReg, s=lambda))[-1]
  } else {
    stop("Your choice is not admisible: choose either OLS or ridge")
  }
  
  # Initialize container for best lambdas
  # and corresponding CV-errors
  lambdaError <- matrix(NA, length(gammaSeq), 2)
  
  # Loop on gamma values
  for (i in 1:length(gammaSeq)) {
    adaGamma <- gammaSeq[i]
    weights <- 1 / abs(initCoefs)^adaGamma
    adaLassoModel <- glmnet::cv.glmnet(X, y, alpha = 1, 
                                       standardize = standardize,
                                       penalty.factor=weights)
    if (lambdaOption == "min") {
      lambdaOpt <- adaLassoModel$lambda.min
    } else if (lambdaOption == "1se") {
      lambdaOpt <- adaLassoModel$lambda.1se
    } else {
      stop("Error: your choice of lambda option should be either min or 1se")
    }
    lambdaError[i, 1] <- lambdaOpt
    lambdaError[i, 2] <- adaLassoModel$cvm[adaLassoModel$lambda == lambdaOpt]
  }
  
  # Find position of minimum CV-error
  posMin <- which.min(lambdaError[, 2])
  
  # Extract best (lambda, gamma)-pair
  bestLambda <- lambdaError[posMin, 1]
  bestGamma <- gammaSeq[posMin]
  minError <- lambdaError[posMin, 2]
  weights <- 1 / abs(initCoefs)^bestGamma
  
  # Get Adaptive Lasso model with best combination
  bestModel <- glmnet::glmnet(X, y, alpha=1,
                              penalty.factor=weights,
                              standardize = standardize)
 
  adaList <- list(AdaModel = bestModel,
                  lambdaOpt = bestLambda,
                  gammaOpt = bestGamma,
                  minError = minError)
  
  # Return model and best (lambda, gamma)-combination 
  return(adaList)
}
