#' Replicate the results of the model selection consistency study
#' 
#' Replicate the results of the Monte Carlo Simulation for the model selection
#' consistency study. The default values are set such that the results you obtain
#' are the same as in our term paper.
#' @param path The path to which you want the graphics to be exported to. 
#' @param scenario Big or small coefficients scenarios to be simulated. Default are both.
#' @param nObs Number of observations.
#' @param p Number of variables.
#' @param nSim Number of simulations.
#' @param nDesigns Number of design matrices.
#' @param sigma Standard deviation of the noise for the regression model.
#' @param adaGamma Gamma parameter for the adaptive Lasso. Default is 2.
#' @param seed Seed for the simulations
#' @return A (Two) pdf-file(s) exported to your path. Moreover, the 
#' results are also saved in a list if you assign a variable 
#' to the function as in the example below.
#' @keywords LassoGroupProject
#' @export
#' @examples results <- ModelSelSim() # It takes around 1 hour
ModelSelSim <- function(path = getwd(), scenario = c("big", "small"), 
                     nObs = 1000, p = 50, nSim = 100, nDesigns = 100, 
                     sigma = 0.1, adaGamma = 2, seed = 481) {

  # Call ggplot
  library("ggplot2")

  # Save current directory to return back
  currentDir <- getwd()

  # Create folder structure if it does not exist
  simPath <- file.path(path, "Model_Selection_Simulation")
  dir.create(simPath, showWarnings = FALSE)

  # Set working directory for exporting the graphs
  setwd(simPath)
    
  # Set seed
  set.seed(seed)
  
  # Number of models = 2: Big vs Small coefficients
  nModels <- length(scenario)
  
  # Coefficients
  bigCoef <- c(-5, 4, 10, -7, 8, 4, 9, -8, -3, 8)
  smallCoef <- c(0.1, -0.2, 0.9, 0.3, -0.8, -0.3, 0.7, -0.4, -0.2, 0.9)
  
  # Pre-processing parameter choices
  coefMat <- as.matrix(bigCoef)
  coefMat <- cbind(coefMat, smallCoef)
  attributes(coefMat)$dimnames[[2]] <- NULL # clean up names
  
  # Initialize containers
  IRmat <- matrix(NA, 2, nDesigns) # Irrepresentable constants
  lassoProbModel <- matrix(NA, 2, nDesigns) # Probability for Lasso selecting the true model
  adaProbModel <- matrix(NA, 2, nDesigns) # Probability for AdaLasso selecting the true model
  stnr <- matrix(NA, 2, nDesigns) # Signal to noise ratio defined
  avgsnr <- numeric(2) # Average signal to noise ratio
  
  # Sample nDesigns Wishart(p, p) covariance matrices
  Xarray <- stats::rWishart(nDesigns, p, diag(rep(1, p)))

  # Loop over models (Big vs Small coefficients)
  tictoc::tic.clearlog() 
  tictoc::tic() # Count time
  for (i in 1:nModels) {
    
    if (scenario[i] == "big" && i == 2) next 
    
    if (scenario[i] == "small" && i == 1) {
      print(paste0("Scenario number: ", 1, " of ", nModels))
      next
    }
    
    print(paste0("Scenario number: ", i, " of ", nModels))
    
    # Vector of non-zero coefficients
    nzCoef <- coefMat[, i]
    
    # Vector of all coefficients
    coefVec <- c(nzCoef, rep(0, (p - 10)))
    
    # Loop over design matrices
    for (j in 1:nDesigns) {
      
      print(paste0("Design matrix number: ", j, " of ", nDesigns))
      
      # Get the Wishart(p, p) covariance matrix sampled before
      covMat <- Xarray[, , j]
      
      # Design matrix: Centered Regressors
      X <- MASS::mvrnorm(n=nObs, mu=rep(0, p), Sigma=covMat)
      
      # Standardized variables to have 
      # sample variance equal to 1
      # Remark: This makes the coefficients
      # comparable for variables in different scales
      X <- t(t(X) / sqrt(colMeans(X^2) - (colMeans(X))^2))
      
      # Get submatrices
      GramS <- t(X[, 1:10]) %*% X[, 1:10] # Gram matrix of the support
      GramSc <- t(X[, 11:p]) %*% X[, 1:10] # Gram matrix of the regressors not in the support
      
      # Calculate Irrepresentable Constant
      IRmat[i, j] <- max(abs(GramSc %*% solve(GramS) %*% sign(nzCoef)))
      
      # True function
      truth <- X %*% coefVec
      
      # Signal-to-noise-Ratio
      stnr[i, j] <- sqrt(sum(truth^2)) / (sqrt(nObs) * sigma)
      
      # True sign function
      trueSign <- sign(coefVec)
      
      # Initialize containers
      lassoModelCheck <- numeric(nSim)
      adaModelCheck <- numeric(nSim)
      
      # Loop over simulations
      for (k in 1:nSim) {
        
        # Generate noise
        noise <- rnorm(nObs, mean=0, sd=sigma)
        
        # Generate data for the response
        y <- truth + noise			
        
        # Standard Lasso
        # LARS algorithm 
        model <- lars::lars(X, y, type = "lasso", normalize = FALSE)
        
        # Coefficient matrix
        estMat <- t(as.matrix(coef(model)))
        signEstMat <- sign(estMat) # Signs of the coefficients
        
        # Check whether there is a model with the correct signs
        lassoModelCheck[k] <- ifelse(
          min(
            colSums(
              ifelse(signEstMat == trueSign, 0, 1))) == 0,
          1, 0)
        
        # Data frame for the regression
        df <- data.frame(y = y, X)
        
        # OLS as a first step for the Adaptive Lasso weights
        olsModel <- lm(y ~ ., data = df)
        
        # Get here OLS Coefficients
        olsCoef <- as.vector(coef(olsModel))[-1]
        
        # Algorithm from Zou (2006): Section 3.5 Computations
        # The LARS algorithm for the adaptive Lasso
        
        # 1. Define x** = x_j / w_j = x_j * olsCoef^2
        Xdoublestar <- t(t(X) * abs(olsCoef)^adaGamma)
        
        # 2. Solve the Lasso problem
        # Adaptive Lasso via the LARS algorithm
        ada.lasso <- lars::lars(Xdoublestar, y, type = "lasso", normalize = FALSE)
        
        # Coefficient matrix
        estMat <- t(as.matrix(coef(ada.lasso)))
        
        # 3. Output = ada.beta = beta / weight = beta * olsCoef^2
        estMat <- estMat * abs(olsCoef)^adaGamma
        signEstMat <- sign(estMat) # Signs of the coefficients
        
        # Check whether there is a model with the correct signs
        adaModelCheck[k] <- ifelse(
          min(
            colSums(
              ifelse(signEstMat == trueSign, 0, 1))) == 0,
          1, 0)
        
      }
      
      # Probabilities for selecting the true model
      lassoProbModel[i, j] <- mean(lassoModelCheck)
      adaProbModel[i, j] <- mean(adaModelCheck)
      
    }
    
    # Average Signal to Noise Ratio
    avgsnr[i] <- mean(stnr[i, ])
    
  }
  tictoc::toc(log=TRUE)

  if ("big" %in% scenario) {
    # Plots
    # Big coefficients: Lasso vs AdaLasso
    IRframe <- data.frame(IR = IRmat[1,], value = lassoProbModel[1,],
                          model = "Lasso")
    IRframe2 <- data.frame(IR = IRmat[1,], value = adaProbModel[1,],
                           model = "Adaptive Lasso")
    IRframe <- rbind(IRframe, IRframe2)

    pt <- ggplot(data = IRframe, aes(x=IR, y=value)) +
      geom_point() +
      labs(x = "Irrepresentable Constant",
           y = "Percentage of correct selected models",
           title = "Lasso vs Adaptive Lasso: Model selection",
           subtitle = "With big coefficients",
           caption = paste0("Average Signal-to-Noise Ratio: ", round(avgsnr[1], 1))) +
      facet_grid(. ~ model) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(filename = "Lasso_AdaLasso_IR_BigCoef.pdf", 
           plot = pt, 
           device = pdf,
           height=9,
           width=15,
           units="cm")
    
    # Save results in a list to return
    ModelSelResults <- list(BigCoef = IRframe)
  }
  
  if ("small" %in% scenario) {
    # Small coefficients: Lasso vs AdaLasso
    IRframe <- data.frame(IR = IRmat[2,], value = lassoProbModel[2,],
                          model = "Lasso")
    IRframe2 <- data.frame(IR = IRmat[2,], value = adaProbModel[2,],
                           model = "Adaptive Lasso")
    IRframe <- rbind(IRframe, IRframe2)
    pt <- ggplot(data = IRframe, aes(x=IR, y=value)) +
      geom_point() +
      labs(x = "Irrepresentable Constant",
           y = "Percentage of correct selected models",
           title = "Lasso vs Adaptive Lasso: Model selection",
           subtitle = "With small coefficients",
           caption = paste0("Average Signal-to-Noise Ratio: ", round(avgsnr[2], 1))) +
      facet_grid(. ~ model) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(filename="Lasso_AdaLasso_IR_SmallCoef.pdf", 
           plot = pt, 
           device = pdf,
           height=9,
           width=15,
           units="cm")
    
    # Save results in a list to return
    if ("big" %in% scenario) {
      ModelSelResults[["SmallCoef"]] <- IRframe
    } else {
      ModelSelResults <- list(SmallCoef = IRframe)
    }
  }
  
  # Save results to disk
  save(ModelSelResults, file="ModelSelResults.RData")
  
  # Return results
  return(ModelSelResults)
  
  # Return to the directory at the beginning
  setwd(currentDir)
}
