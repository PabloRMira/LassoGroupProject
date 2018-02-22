#' Replicate the results of the prediction contest simulation study
#' 
#' Replicate the results of the Monte Carlo Simulation for the prediction
#' contest. The default values are set such that the results you obtain
#' are the same as in our term paper.
#' @param path The path to which you want the results to be exported to. 
#' @param pList Vector with the choices of number of variables.
#' @param sigmaList Vector with the choices of noise levels.
#' @param nObs Number of observations.
#' @param nSim Number of simulations.
#' @param adaGamma Gamma parameter for the adaptive Lasso. Default is 2.
#' @param seed Seed for the simulations
#' @return (Number of variable choices x Number of noise levels) plots and 
#' tables in LaTeX and html.
#' @details This function creates a new directory in the directory specified
#' in the 'path' argument named 'Prediction_Contest'. This directory in turn
#' include two subdirectories, 'Plots' and 'Tables'. Finally, in the 'Tables'
#' directory are included two directories more, 'HTML' for the result tables
#' in html-format and 'LaTeX' for the result tables in LaTeX format.
#' @keywords LassoGroupProject
#' @export
PredContest <- function(path = getwd(), pList = c(50, 100),
                        sigmaList = c(1, 5, 25), nObs = 1000,
                        nSim = 100, adaGamma = 2, seed = 512) {

  # Call ggplot
  library("ggplot2")

  # Create folder structure if it does not exist
  simPath <- file.path(path, "Prediction_Contest")
  plotsPath <- file.path(simPath, "Plots")
  tablesPath <- file.path(simPath, "Tables")
  htmlPath <- file.path(tablesPath, "HTML")
  latexPath <- file.path(tablesPath, "LaTeX")
  
  dir.create(simPath, showWarnings = FALSE)
  dir.create(plotsPath, showWarnings = FALSE)
  dir.create(tablesPath, showWarnings = FALSE)
  dir.create(htmlPath, showWarnings = FALSE)
  dir.create(latexPath, showWarnings = FALSE)

  # Set seed
  set.seed(seed)
  
  # Number of variables scenarios: p = 50 vs p = 100
  nModels <- length(pList)
  
  # Number of noise scenarios: sigma in {1, 3, 25}
  nSigmas <- length(sigmaList)

  # Non-zero coefficients  
  nzCoef <- c(-5, 4, 10, -7, 8, 4, 9, -8, -3, 8)
  
  # Loop over scenarios (low vs high sparseness)
  tictoc::tic.clearlog() 
  tictoc::tic() # Count time
  for(i in 1:nModels) {
    
    print(paste0("Variable number scenario: ", i, " of ", nModels))
    
    # Number of variables
    p <- pList[i]
    
    # Complete coefficient vector
    coefVec <- c(nzCoef, rep(0, p-10))
    
    # Sampled covariance matrix from Wishart distribution
    covMat <- stats::rWishart(1, p, diag(rep(1, p)))[,,1]
    
    # Design matrix: Centered Regressors
    X <- MASS::mvrnorm(n=nObs, mu=rep(0, p), Sigma=covMat)
    
    # Standardized variables to have 
    # sample variance equal to 1
    # Remark: This makes the coefficients
    # comparable for variables in different scales
    X <- t(t(X) / sqrt(colMeans(X^2) - (colMeans(X))^2))
    
    # Initialize containers
    # Fit container
    methodNumber <- 5 # Number of methods to evaluate
    modelFit <- array(NA, c(nObs, nSim, methodNumber))
    
    # Results as array
    resultsArray <- array(NA, c(methodNumber, 3))
    
    # True function
    truth <- X %*% coefVec
    
    # Get submatrices
    GramS <- t(X[, 1:10]) %*% X[, 1:10] # Gram matrix of the support
    GramSc <- t(X[, 11:p]) %*% X[, 1:10] # Gram matrix of the regressors not in the support
    
    # Irrepresentable constant
    IR <- max(abs(GramSc %*% solve(GramS) %*% sign(nzCoef)))
    
    # Loop over noise scenarios
    for (j in 1:nSigmas) {
      
      print(paste0("Sigma scenario: ", j, " of ", nSigmas))
      sigma <- sigmaList[j]
      
      # Signal-to-noise Ratio
      # definition in Bühlmann and van der Geer
      snratio <- sqrt(sum(truth^2)) / (sqrt(nObs) * sigma)
      
      # Loop over simulations
      for (k in 1:nSim) {
        
        if (k %% 10 == 0) {
          print(paste0("Simulation number: ", k, " of ", nSim))
        }
        
        # Generate noise
        noise <- rnorm(nObs, mean=0, sd=sigma)
        
        # Generate data for the response
        y <- truth + noise			
        
        # Data frame for the regression
        df <- data.frame(y = y, X)
        
        ############################################
        ############### OLS ########################
        ############################################
        
        # Fitting the model    
        olsModel <- lm(y ~ ., data = df)
        
        # Model fit
        modelFit[, k, 1] <- olsModel$fitted.values
        
        # Get here OLS Coefficients for AdaLasso
        olsCoef <- as.vector(coef(olsModel))[-1]
        
        #######################################
        ############## LASSO ##################
        #######################################
        
        # Standard Lasso with Coordinate Descent (glmnet) 
        lassoModel <- glmnet::cv.glmnet(X, y, alpha = 1, 
                                        standardize = FALSE)
        lambdaOpt <- lassoModel$lambda.min
        
        # Lasso fit
        modelFit[, k, 2] <- predict(lassoModel, 
                                    newx=X, s=lambdaOpt)
        
        #############################################
        ############ Adaptive Lasso #################
        #############################################
        
        weights <- 1 / olsCoef^adaGamma
        adaLassoModel <- glmnet::cv.glmnet(X, y, alpha = 1, 
                                           standardize = FALSE,
                                           penalty.factor=weights)
        lambdaOpt <- adaLassoModel$lambda.min
        
        # AdaLasso Fit
        modelFit[, k, 3] <- predict(adaLassoModel, 
                                    newx=X, s=lambdaOpt)
        
        ##########################################################
        ############## Forward Stepwise Regression ###############
        ##########################################################
        
        null <- lm(y ~ 1, data = df)
        full <- lm(y ~ ., data = df)
        
        forwardModel <- stats::step(null, scope=list(lower=null, upper=full), 
                                    direction="forward",
                                    trace = 0)
        modelFit[, k, 4] <- forwardModel$fitted.values
        
        ##########################################################
        ############# Backward Stepwise Regression ###############
        ##########################################################
        
        backwardModel <- stats::step(full, data = df, direction="backward",
                                     trace  = 0)
        modelFit[, k, 5] <- backwardModel$fitted.value
        
      }
      
      # Save results
      for (s in 1:5) {
        resultsArray[s, 1] <- LassoGroupProject::ABias_fun(modelFit[,,s], truth)
        resultsArray[s, 2] <- LassoGroupProject::AVar_fun(modelFit[,,s])
        resultsArray[s, 3] <- LassoGroupProject::AMSPE_fun(modelFit[,,s], truth)
      }
      
      # Table: Sparse-Model
      sparseTable <- as.data.frame(resultsArray)
      names(sparseTable) <- c("Squared ABias", "AVariance", "AMSPE")
      rownames(sparseTable) <- c("OLS", "Lasso", "Adaptive Lasso", "Forward Stepwise",
                                 "Backward Stepwise")
      
      # Save in corresponding folder
      setwd(htmlPath)
      
      # Output HTML
      stargazer::stargazer(sparseTable, 
                           type="html", 
                           out=paste0("PredContest_Sigma_", sigma, "_p_", p, ".html"), 
                           title="Prediction Contest", 
                           summary=FALSE,
                           notes=c(paste0("Number of Observations: ", nObs),
                                   paste0("Number of Variables: ", p),
                                   paste0("Number of non-zero coefficients: ", 10),
                                   paste0("Number of Simulations: ", nSim),
                                   paste0("Sigma: ", sigma),
                                   paste0("Signal-to-Noise Ratio: ", round(snratio[1], 1)),
                                   paste0("Irrepresentable Constant: ", round(IR, 1))))
      
      setwd(latexPath)
      
      # Output LaTeX
      stargazer::stargazer(sparseTable, 
                           type="latex", 
                           out=paste0("PredContest_Sigma_", sigma, "_p_", p, ".html"), 
                           title="Prediction contest", 
                           summary=FALSE,
                           notes=c(paste0("Number of Observations: ", nObs),
                                   paste0("Number of Variables: ", p),
                                   paste0("Number of non-zero coefficients: ", 10),
                                   paste0("Number of Simulations: ", nSim),
                                   paste0("Sigma: ", sigma),
                                   paste0("Signal-to-Noise Ratio: ", round(snratio, 1)),
                                   paste0("Irrepresentable Constant: ", round(IR, 1))))
      
      setwd(plotsPath)
      
      # Plot results
      df <- sparseTable[-3]
      df$model <- row.names(df)
      mdf <- reshape::melt(df)
      mdf$modelOrd <- factor(mdf$model, levels=c("OLS", "Lasso", "Adaptive Lasso",
                                                 "Forward Stepwise", "Backward Stepwise"))
      
      pt <- ggplot(data = mdf, aes(x=modelOrd, y=value, fill=variable)) +
        geom_col() +
        labs(title="AMSPE Decomposition",
             subtitle=paste0("Number of Variables: ", p, ", SNR: ", round(snratio, 1)),
             caption=paste0("Irrepresentable Constant: ", round(IR, 2)),
             x="",
             y="") +
        theme(plot.title=element_text(hjust=0.5),
              plot.subtitle=element_text(hjust=0.5),
              legend.position="bottom",
              legend.title.align=0.5) +
        guides(fill=guide_legend(title="Components",
                                 title.position="top"))
      
      ggsave(filename=paste0("PredContest_Sigma_", sigma, "_p_", p, ".pdf"),
             plot=pt,
             scale=1.3)
    }
  }
  tictoc::toc(log=TRUE)
}
