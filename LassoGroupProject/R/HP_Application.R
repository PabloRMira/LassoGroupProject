#' Replicate the results of the house prices application
#' 
#' Replicate the results of the real data application with the 
#' house prices data set.
#' @param path The path to which you want the results to be exported to. 
#' @param adaGamma Gamma parameter for the adaptive Lasso. Default is 2.
#' @param mccv Should Monte Carlo Cross Validation be applied for evaluation? Default is TRUE.
#' @param cvIter Number of iterations for the Monte Carlo Cross Validation. Default is 100.
#' @param seed Seed for the simulations.
#' @return Result tables and plots.
#' @details This function creates a new directory, 'Application_HP' where
#' the results will be saved.
#' @keywords LassoGroupProject
#' @export
HP_Application <- function(path = getwd(), 
                           adaGamma = 2,
                           mccv = TRUE,
                           cvIter = 100,
                           seed = 123) {

  # Call ggplot
  library("ggplot2", suppressPackageStartupMessages())
  
  # Create folder structure if it does not exist
  appPath <- file.path(path, "Application_HP")
  dir.create(appPath, showWarnings = FALSE)

  # Set seed
  set.seed(seed)
  
  # Define house sale price as dependent variable
  y <- as.vector(as.matrix(hp_data[32]))
  
  # Define all other varaibles as explanatory variables
  X <- as.matrix(hp_data[-32])

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
  model <- glmnet::cv.glmnet(X, y, family="gaussian",
                                  alpha=0,
                             standardize=TRUE)
  lambda <- model$lambda.1se
  
  # Extract ridge coefficients
  ridgeCoefs <- as.vector(coef(model, s=lambda))[-1]
  
  # Use ridge coefficients to compute the weights
  # for the adaptive Lasso
  weights <- 1 / abs(ridgeCoefs)^adaGamma
 
  # Cross-Validation
  model <- glmnet::cv.glmnet(X, y, family="gaussian",
                             alpha=1,
                             penalty.factor=weights,
                             standardize=TRUE)
  lambda <- model$lambda.1se

  # Extract cross-validated mean squared
  # prediction error (MSPE)
  adaMSPE <- model$cvm[model$lambda == lambda]

  # Extract coefficients
  coeffi <- coef(model, s=lambda)
  adaCoeffi <- as.vector(coeffi)
  
  # Get normalized coefficients by multiplying 
  # each coefficient with the standard deviation
  # of the explanatory variables
  adaCoeffiNorm <- c(adaCoeffi[1],
                     adaCoeffi[-1] * sqrt(colMeans(X^2) - (colMeans(X))^2))

  adaSelVarPos <- coeffi@i + 1
  adaSelVarVal <- adaCoeffiNorm[adaSelVarPos][-1]
  adaSelVarNum <- length(adaSelVarVal)
  adaSelVarName <- coeffi@Dimnames[[1]][adaSelVarPos][-1]

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
      model <- glmnet::cv.glmnet(X, y, family="gaussian",
                                 alpha=0,
                                 standardize=TRUE)
      lambda <- model$lambda.1se
      
      # Extract ridge coefficients
      ridgeCoefs <- as.vector(coef(model, s=lambda))[-1]
      
      # Use ridge coefficients to compute the weights
      # for the adaptive Lasso
      weights <- 1 / abs(ridgeCoefs)^adaGamma
      
      # Cross-Validation
      model <- glmnet::cv.glmnet(X, y, family="gaussian",
                                 alpha=1,
                                 penalty.factor=weights,
                                 standardize=TRUE)
      lambda <- model$lambda.1se
      adaCVASPE[i] <- model$cvm[model$lambda == lambda]
    }
    # Kernel density plot from Monte Carlo Cross Validation
    dfMCCV <- data.frame(Lasso = lassoCVASPE, AdaLasso = adaCVASPE)
    names(dfMCCV) <- c("Lasso", "Adaptive Lasso")
    mdfMCCV <- reshape::melt(dfMCCV)
    names(mdfMCCV) <- c("model" , "value")
    pt <- ggplot(data = mdfMCCV, aes(x = value, colour=model, fill=model)) +
      geom_density(alpha=0.7) +
      labs(title="House Prices Application: Monte Carlo Cross Validation",
           x="MCCVASPE",
           y="Kernel Density Estimator") +
      theme(plot.title=element_text(hjust=0.5)) +
      guides(fill = guide_legend(title=""),
             colour = guide_legend(title=""))
    
    setwd(appPath)
    ggsave(filename="HP_Application_MCCV.pdf",
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

  setwd(appPath)
  if (mccv == TRUE) {
    # Monte Carlo Cross-Validation
    lassoMCCVASPE <- mean(lassoCVASPE)
    adaMCCVASPE <- mean(adaCVASPE)
    # Results
    stargazer::stargazer(fullTable,
                         type="html",
                         out="HP_Results_MCCV.html",
                         title="House Prices Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso MCCVASPE: ", round(lassoMCCVASPE, 4)),
                                 paste0("Adaptive MCCVASPE: ", round(adaMCCVASPE, 4))))
    
    stargazer::stargazer(fullTable,
                         type="latex",
                         out="SP_Results.tex",
                         title="S&P Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso MCCVASPE: ", round(lassoMCCVASPE, 4)),
                                 paste0("Adaptive MCCVASPE: ", round(adaMCCVASPE, 4))))
  } else {
    # Plain vanilla Cross-Validation
    stargazer::stargazer(fullTable,
                         type="html",
                         out="HP_Results.html",
                         title="House Prices Application",
                         summary=FALSE,
                         digits=3,
                         notes=c(paste0("Lasso CVASPE: ", round(lassoMSPE, 4)),
                                 paste0("Adaptive Lasso CVASPE: ", round(adaMSPE, 4))))
    
    stargazer::stargazer(fullTable,
                         type="latex",
                         out="HP_Results.tex",
                         title="House Prices Application",
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

  ggsave(filename="HP_Application_Coefficients.pdf",
         plot=pt,
         height=9,
         width=15,
         units="cm")
}
