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
#' @param adaList Sequence of gamma parameters for the 
#' adaptive Lasso to cross-validate.
#' @param seed Seed for the simulations.
#' @return (Number of variable choices x Number of noise levels) plots and 
#' tables in LaTeX and html. Moreover, the results are also saved in a list
#' if you assign a variable to the function as in the example below. The 
#' results are also exported to the directory specified in the 'path' argument.
#' @details This function creates a new directory in the directory specified
#' in the 'path' argument named 'Prediction_Contest'. This directory in turn
#' include two subdirectories, 'Plots' and 'Tables'. Finally, in the 'Tables'
#' directory are included two directories more, 'HTML' for the result tables
#' in html-format and 'LaTeX' for the result tables in LaTeX format.
#' @keywords LassoGroupProject
#' @export
#' @examples results <- PredContest() # It takes long time with default arguments (around 5 hours)
PredContest <- function(path = getwd(), 
                        pList = c(50, 100),
                        sigmaList = c(1, 5, 25), 
                        nObs = 500,
                        nSim = 100, 
                        adaList = c(0.01, 0.1, 0.5, 1, 2, 5, 7), 
                        seed = 823) {

  # Call ggplot
  library("ggplot2")

  # Save current directory to return back
  currentDir <- getwd()
  
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

  # Number of methods
  methodNumber <- 5
  
  # Initialize container list for the results
  PredResults <- list()
  
  # Initialize counter for the list index
  countList <- 0
  
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
      
      # Update counter for list
      countList <- countList + 1
      
      # Signal-to-noise Ratio
      # definition in Bühlmann and van der Geer
      snratio <- sqrt(sum(truth^2)) / (sqrt(nObs) * sigma)

      # Initialize container for the time
      methodTime <- matrix(NA, nSim, methodNumber)

      # Initialize containers
      # Fit container
      modelFit <- array(NA, c(nObs, nSim, methodNumber))
      
      # Results as array
      resultsArray <- matrix(NA, methodNumber, 4)

      # Initialize container for optimal gamma sequence
      gammaVec <- numeric(nSim)
            
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
        startTimer <- Sys.time()
        olsModel <- lm(y ~ ., data = df)
        endTimer <- Sys.time()
        
        methodTime[k, 1] <- round(as.vector(endTimer - startTimer) / 60, 3)
        
        # Model fit
        modelFit[, k, 1] <- olsModel$fitted.values
        
        #######################################
        ############## LASSO ##################
        #######################################
        
        # Standard Lasso with Coordinate Descent (glmnet) 
        startTimer <- Sys.time()
        lassoModel <- glmnet::cv.glmnet(X, y, alpha = 1, 
                                        standardize = FALSE)
        endTimer <- Sys.time()

        methodTime[k, 2] <- round(as.vector(endTimer - startTimer) / 60, 3)
                
        lambdaOpt <- lassoModel$lambda.min
        
        # Lasso fit
        modelFit[, k, 2] <- predict(lassoModel, 
                                    newx=X, s=lambdaOpt)
        
        #############################################
        ############ Adaptive Lasso #################
        #############################################

        startTimer <- Sys.time()
        adaLassoModel <- LassoGroupProject::AdalassoCV(X, y, 
                                                       gammaSeq = adaList,
                                                       standardize=FALSE)
        endTimer <- Sys.time()
        
        methodTime[k, 3] <- round(as.vector(endTimer - startTimer) / 60, 3)
        
        lambdaOpt <- adaLassoModel$lambdaOpt
        adaModel <- adaLassoModel$AdaModel
        gammaVec[k] <- adaLassoModel$gammaOpt
        
        # AdaLasso Fit
        modelFit[, k, 3] <- predict(adaModel, 
                                    newx=X, s=lambdaOpt)
        
        ##########################################################
        ############## Forward Stepwise Regression ###############
        ##########################################################
        
        startTimer <- Sys.time()
        null <- lm(y ~ 1, data = df)
        full <- lm(y ~ ., data = df)
        forwardModel <- stats::step(null, scope=list(lower=null, upper=full), 
                                    direction="forward",
                                    trace = 0)
        endTimer <- Sys.time()
        
        methodTime[k, 4] <- round(as.vector(endTimer - startTimer) / 60, 3)
        
        modelFit[, k, 4] <- forwardModel$fitted.values
        
        ##########################################################
        ############# Backward Stepwise Regression ###############
        ##########################################################
        
        startTimer <- Sys.time()
        backwardModel <- stats::step(full, data = df, direction="backward",
                                     trace  = 0)
        endTimer <- Sys.time()
        
        methodTime[k, 5] <- round(as.vector(endTimer - startTimer) / 60, 3)

        modelFit[, k, 5] <- backwardModel$fitted.value
        
      }
      
      # Save results
      for (s in 1:5) {
        resultsArray[s, 1] <- LassoGroupProject::ABias_fun(modelFit[,,s], truth)
        resultsArray[s, 2] <- LassoGroupProject::AVar_fun(modelFit[,,s])
        resultsArray[s, 3] <- LassoGroupProject::AMSPE_fun(modelFit[,,s], truth)
        resultsArray[s, 4] <- mean(methodTime[, s])
      }
      
      # Average gamma
      avgGamma <- round(mean(gammaVec), 3)
      
      # Table: Sparse-Model
      sparseTable <- as.data.frame(resultsArray)
      names(sparseTable) <- c("Squared ABias", "AVariance", "AMSPE", "Time consumed")
      rownames(sparseTable) <- c("OLS", "Lasso", "Adaptive Lasso", "Forward Stepwise",
                                 "Backward Stepwise")
      
      # Save results in a list
      PredResults[[countList]] <- sparseTable
      names(PredResults)[countList] <- paste0("Sigma_", sigma, "_p_", p)
      
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
                                   paste0("Average optimal gamma: ", avgGamma),
                                   paste0("Irrepresentable Constant: ", round(IR, 1))))
      
      setwd(latexPath)
      
      # Output LaTeX
      stargazer::stargazer(sparseTable, 
                           type="latex", 
                           out=paste0("PredContest_Sigma_", sigma, "_p_", p, ".tex"), 
                           title="Prediction contest", 
                           summary=FALSE,
                           notes=c(paste0("Number of Observations: ", nObs),
                                   paste0("Number of Variables: ", p),
                                   paste0("Number of non-zero coefficients: ", 10),
                                   paste0("Number of Simulations: ", nSim),
                                   paste0("Sigma: ", sigma),
                                   paste0("Signal-to-Noise Ratio: ", round(snratio, 1)),
                                   paste0("Average optimal gamma: ", avgGamma),
                                   paste0("Irrepresentable Constant: ", round(IR, 1))))
      
      setwd(plotsPath)
      
      # Plot results
      df <- sparseTable[c(-3, -4)]
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
             height=9,
             width=17,
             units="cm")
    }
  }
  tictoc::toc(log=TRUE)
  
  # Save results to disk
  setwd(simPath)
  save(PredResults, file="PredResults.RData")
  
  # Return results
  return(PredResults)
  
  # Return back to initial directory
  setwd(currentDir)
}
