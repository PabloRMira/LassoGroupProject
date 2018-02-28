#' Replicate the results of the stock prices application
#' 
#' Replicate the results of the real data application with the 
#' S&P data set.
#' @param path The path to which you want the results to be exported to. 
#' @param adaList Gamma parameter for the adaptive Lasso. Default is 1.
#' @param mccv Should Monte Carlo Cross Validation be applied for evaluation? Default is TRUE.
#' @param cvIter Number of iterations for the Monte Carlo Cross Validation. Default is 100.
#' @param seed Seed for the simulations.
#' @return Result tables and plots.
#' @details This function creates a new directory, 'Application_SP' where
#' the results will be saved.
#' @keywords LassoGroupProject
#' @export
#' @examples SP_results <- SP_Application()
SP_Application <- function(path = getwd(),
                           adaList = c(0.01, 0.1, 0.5, 1, 2, 3),
                           mccv = TRUE,
                           cvIter = 100,
                           seed = 123) {
  
  # Call ggplot
  library("ggplot2", suppressPackageStartupMessages())

  # Get current directory to return back
  currentDir <- getwd()
    
  # Create folder structure if it does not exist
  appPath <- file.path(path, "Application_SP")
  dir.create(appPath, showWarnings = FALSE)
  
  # Set seed
  set.seed(seed)
  
  # Define stock price as dependent variable
  y <- as.vector(as.matrix(sp_data$avg_price))
  
  # Define all other varaibles as explanatory variables
  X <- as.matrix(sp_data[, !(names(sp_data) == "avg_price")])
  
  #############
  ### Lasso ###
  #############
  
  # Cross-Validation
  model <- glmnet::cv.glmnet(X, y, family="gaussian", 
                             alpha=1, standardize=TRUE)
  lambda <- model$lambda.1se
  
  # Extract cross-validated mean squared
  # prediction error (MSPE)
  lassoMSPE <- model$cvm[model$lambda == lambda]
  
  # Extract coefficients
  coeffi <- coef(model, s=lambda)
  lassoCoeffi <- as.vector(coeffi)
  
  # Get normalized coefficients by multiplying 
  # each coefficient with the standard deviation
  # of the explanatory variables
  lassoCoeffiNorm <- c(lassoCoeffi[1],
                       lassoCoeffi[-1] * sqrt(colMeans(X^2) - (colMeans(X))^2))
  
  lassoSelVarPos <- coeffi@i + 1
  lassoSelVarVal <- lassoCoeffiNorm[lassoSelVarPos][-1]
  lassoSelVarNum <- length(lassoSelVarVal)
  lassoSelVarName <- coeffi@Dimnames[[1]][lassoSelVarPos][-1]
  
  #######################
  ### Adaptive Lasso ###
  #######################
  
  set.seed(seed)
  
  # Use ridge regression as first stage estimator
  # for the adaptive Lasso since there is 
  # multicollinearity in the data set and hence
  # OLS breaks down in such an environment
  # Cross-Validation
  adaLassoModel <- LassoGroupProject::AdalassoCV(X, y, 
                                                 gammaSeq = adaList,
                                                 standardize=TRUE,
                                                 lambdaOption = "1se",
                                                 initEst = "ridge")
  lambda <- adaLassoModel$lambdaOpt
  model <- adaLassoModel$AdaModel

  # Extract cross-validated mean squared
  # prediction error (MSPE)
  adaMSPE <- adaLassoModel$minError
  
  # Extract coefficients
  coeffi <- coef(model, s=lambda)
  adaCoeffi <- as.vector(coef(model, s=lambda))
  
  # Get normalized coefficients by multiplying 
  # each coefficient with the standard deviation
  # of the explanatory variables
  adaCoeffiNorm <- c(adaCoeffi[1],
                     adaCoeffi[-1] * sqrt(colMeans(X^2) - (colMeans(X))^2))
  
  adaSelVarPos <- coeffi@i + 1
  adaSelVarVal <- adaCoeffiNorm[adaSelVarPos][-1]
  adaSelVarNum <- length(adaSelVarVal)
  adaSelVarName <- coeffi@Dimnames[[1]][adaSelVarPos][-1]

  setwd(appPath)
    
  if (mccv == TRUE) {
    # Monte Carlo Cross Validation for Lasso and adaptive Lasso
    # Evaluation criterion: 
    # Cross-Validated Average Squared Prediction Error (CVASPE)
    set.seed(seed)
    lassoCVASPE <- numeric(cvIter)
    adaCVASPE <- numeric(cvIter)
    for (i in 1:cvIter) {
      # Print information
      print(paste0("Iteration number: ", i, " of ", cvIter))
      
      # Lasso
      model <- glmnet::cv.glmnet(X, y, family="gaussian", 
                                 alpha=1, standardize=TRUE)
      lambda <- model$lambda.1se
      lassoCVASPE[i] <- model$cvm[model$lambda == lambda]
      
      # Adaptive Lasso
      # First stage: Ridge regression
      # Cross-Validation
      adaLassoModel <- LassoGroupProject::AdalassoCV(X, y, 
                                                     gammaSeq = adaList,
                                                     standardize=TRUE,
                                                     lambdaOption = "1se",
                                                     initEst = "ridge")
      lambda <- adaLassoModel$lambdaOpt
      model <- adaLassoModel$AdaModel
      adaCVASPE[i] <- adaLassoModel$minError
    }
    # Kernel density plot from Monte Carlo Cross Validation
    dfMCCV <- data.frame(Lasso = lassoCVASPE, AdaLasso = adaCVASPE)
    names(dfMCCV) <- c("Lasso", "Adaptive Lasso")
    mdfMCCV <- reshape::melt(dfMCCV)
    names(mdfMCCV) <- c("model" , "value")
    pt <- ggplot(data = mdfMCCV, aes(x = value, colour=model, fill=model)) +
      geom_density(alpha=0.7) +
      labs(title="Stock Price Data Application: Monte Carlo Cross Validation",
           x="ASPE",
           y="Kernel Density Estimator") +
      theme(plot.title=element_text(hjust=0.5)) +
      guides(fill = guide_legend(title=""),
             colour = guide_legend(title=""))
    
    ggsave(filename="SP_Application_MCCV.pdf",
           plot=pt,
           height=9,
           width=15,
           units="cm")
  }
  
  # Coefficient Tables
  lassoTable <- data.frame(name = lassoSelVarName,
                           val = lassoSelVarVal)
  adaTable <- data.frame(name = adaSelVarName,
                         val = adaSelVarVal)
  fullTable <- merge(x = lassoTable, y = adaTable,
                     by = "name", all=TRUE)
  names(fullTable) <- c("Variable name",
                        "Lasso",
                        "Adaptive Lasso")
  # Substitute points by spaces in name
  fullTable[[1]] <- gsub("[.]", " ", fullTable[[1]])

  # Save
  
  if (mccv == TRUE) {
    # Monte Carlo Cross-Validation
    lassoMCCVASPE <- mean(lassoCVASPE)
    adaMCCVASPE <- mean(adaCVASPE)
    # Results
    stargazer::stargazer(fullTable,
                         type="html",
                         out="SP_Results_MCCV.html",
                         title="Stock Price Data Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso MCCVASPE: ", round(lassoMCCVASPE, 4)),
                                 paste0("Adaptive Lasso MCCVASPE: ", round(adaMCCVASPE, 4))))
    
    stargazer::stargazer(fullTable,
                         type="latex",
                         out="SP_Results.tex",
                         title="Stock Price Data Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso MCCVASPE: ", round(lassoMCCVASPE, 4)),
                                 paste0("Adaptive Lasso MCCVASPE: ", round(adaMCCVASPE, 4))))
  } else {
    # Plain vanilla Cross-Validation
    stargazer::stargazer(fullTable,
                         type="html",
                         out="SP_Results.html",
                         title="S&P Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso CVASPE: ", round(lassoMSPE, 4)),
                                 paste0("Adaptive Lasso CVASPE: ", round(adaMSPE, 4))))
    
    stargazer::stargazer(fullTable,
                         type="latex",
                         out="SP_Results.tex",
                         title="Stock Price Data Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso CVASPE: ", round(lassoMSPE, 4)),
                                 paste0("Adaptive Lasso CVASPE: ", round(adaMSPE, 4))))
  }
  
  # Coefficient plots
  
  # Melt data
  mdf <- reshape::melt(fullTable)
  names(mdf) <- c("variable" , "model", "value")
  
  pt <- ggplot(data=mdf, aes(x=variable, y=value)) +
    geom_bar(stat="identity") +
    coord_flip() +
    facet_grid(. ~ model) +
    labs(title="Lasso vs Adaptive Lasso: Non-zero coefficients",
         x="Variable",
         y="Coefficient") +
    theme(plot.title=element_text(hjust=0.5))
  
  ggsave(filename="SP_Application_Coefficients.pdf",
         plot=pt,
         height=9,
         width=15,
         units="cm")
  
  # Save results in a list and to disk
  SP_results <- list(fullTable)
  
  if (mccv==TRUE) {
    SP_results[["MCCV"]] <- data.frame(Lasso = lassoMCCVASPE, 
                                       Adalasso= adaMCCVASPE)
  } else {
    SP_results[["CV"]] <- data.frame(Lasso = lassoMSPE,
                                     Adalasso = adaMSPE)
  }
  
  save(SP_results, file="SP_results.RData")
  
  # Return results
  return(SP_results)
  
  # Return back
  setwd(currentDir)
}
