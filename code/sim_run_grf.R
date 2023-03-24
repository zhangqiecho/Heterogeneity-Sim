#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

## NOTE: We can run this program locally from the command line, and thus leave the package 
## installation as is. Once we run on cluster, need to change library location.

## where is the library?
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","skimr","here","mvtnorm","latex2exp","earth",
   "readxl","VGAM", "ranger","xgboost","mgcv","glmnet","NbClust","factoextra",
  "SuperLearner", "AIPW", "dplyr", "cluster", "ggplot2", "estimatr","grf")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github('susanathey/causalTree')
library(causalTree)

args = c(5,1000)

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

number_sims <- as.numeric(args[1])
cat(paste("Number of Simulations:", number_sims, "\n"))
sample_size <- as.numeric(args[2])
cat(paste("Sample Size for Each Sim:", sample_size), "\n")

## FUNCTION
cluster_sim <- function(nsim, sample_size){
  
  n = sample_size
  p = 5
  set.seed(nsim)
  
  ## CONFOUNDERS
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  c     <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESIGN MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  parms3 <- rep(3,5)
  parms4 <- rep(log(1.5),5)
  beta   <- parms3; beta <- c(120, beta)
  
  # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
  piMatT  <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  piMatT2 <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  theta   <- c(-.5,parms4) # assumed coefficients of propensity score model
  theta2  <- c(-.5,parms4)
  
  # PROPENSITY SCORE MODEL
  ## pi is an actual parameter in base R, so good practice not to overwrite
  pi_x   <- expit(piMatT%*%theta)
  pi_m <- expit(piMatT2%*%theta2)
  x    <- rbinom(n,1,pi_x)
  m    <- rbinom(n,1,pi_m)
  
  # OUTCOME MODEL: EXPOSURE VALUE UNDER M == 0 IS 7.5; VALUE UNDER M = 1 IS -7.5
  y <- 80 + 7.5*x + 7.5*m - 15*x*m + rnorm(length(x),0,10)

  dat <- data.frame(y,m,c,x)
  names(dat) <- c("y","m",paste0("c",1:5),"x")
  
  ## DATA GENERATION FINISHED; DATA ANALYSIS START
  
  # TRUE MODEL
  mod_true <- lm(y ~ x + m + x*m + c1 + c2 + c3 + c4 + c5, data = dat)
  c_matrix1 <- c(0,1,0,0,0,0,0,0,0) # ATE among m = 0
  c_matrix2 <- c(0,1,0,0,0,0,0,0,1) # ATE among m = 1; NB: r puts interaction coefficient at end
  
  ate_m0 <- c(coef(mod_true) %*% c_matrix1, sqrt(t(c_matrix1) %*% vcov(mod_true) %*% c_matrix1))
  
  ate_m1 <- c(coef(mod_true) %*% c_matrix2, sqrt(t(c_matrix2) %*% vcov(mod_true) %*% c_matrix2))
  
  # ## split sample into training and testing
  # train_fraction <- .8  # Use train_fraction % of the dataset to train our models
  # n <- dim(dat)[1]
  # train_idx <- sample.int(n, replace=F, size=floor(n*train_fraction))
  # df_train <- dat[train_idx,]
  # df_test  <- dat[-train_idx,]
  # 
  # ## further split for selection and prediction
  # # Diving the data 40%-40%-20% into splitting, estimation and validation samples
  # split_size <- floor(nrow(df_train) * 0.5)
  # split_idx <- sample(nrow(df_train), replace=FALSE, size=split_size)
  # 
  # # Make the splits
  # df_split <- df_train[split_idx,]
  # df_est <- df_train[-split_idx,]
  # 
  # # regression formula
  # covariate_names <- names(dat)[-1]
  # reg_formula <- paste("y ~", paste(covariate_names, collapse = " + "))
  # 
  # cf <- causal_forest(
  #   X = as.matrix(df_train[,covariate_names]),
  #   Y = df_train$y,
  #   W = df_train$x,
  #   num.trees=1000
  #   )
  
  # Qi: here are your tasks
  ## Ultimate goal: get one mean difference for M = 0, and one for M = 1
  
  W_forest = regression_forest(
    X = as.matrix(data.frame(c,m)), 
    Y = as.matrix(x),
    num.trees = 2000,
    sample.weights = NULL,
    clusters = NULL,
    equalize.cluster.weights = FALSE,
    sample.fraction = 0.5,
    min.node.size = 5,
    honesty = TRUE,
    honesty.fraction = 0.5,
    honesty.prune.leaves = TRUE,
    alpha = 0.05,
    imbalance.penalty = 0
    )
  
  W.hat = predict(W_forest)$predictions
  
  Y_forest = regression_forest(
    X = as.matrix(data.frame(c,m)), 
    Y = as.matrix(y),
    num.trees = 2000,
    sample.weights = NULL,
    clusters = NULL,
    equalize.cluster.weights = FALSE,
    sample.fraction = 0.5,
    min.node.size = 5,
    honesty = TRUE,
    honesty.fraction = 0.5,
    honesty.prune.leaves = TRUE,
    alpha = 0.05,
    imbalance.penalty = 0
    )
  
  Y.hat = predict(Y_forest)$predictions
  
  cf <- causal_forest(
    X = data.frame(c,m), 
    Y = y, 
    W = x,
    Y.hat = Y.hat, 
    W.hat = W.hat, 
    num.trees = 2000, 
    tune.parameters = "all"
    )
  
  # ggplot(data.frame(W.hat = cf$W.hat, W = factor(cf$W.orig))) +
  #   geom_histogram(aes(x = W.hat, y = stat(density), fill = x), alpha=0.3, position = "identity") +
  #   geom_density(aes(x = W.hat, color = x)) +
  #   xlim(0,1) +
  #   labs(title = "Causal forest propensity scores",
  #        caption = "The propensity scores are learned via GRF's regression forest")

  ate_m1_grf <- average_treatment_effect(cf, subset = m == 1) # ate = -6.8 when m = 1
  ate_m0_grf <- average_treatment_effect(cf, subset = m == 0) # ate = 7.7 when m = 0

  res <- data.frame(
    rbind(
      ate_m0,
      ate_m1,
      ate_m0_grf,
      ate_m1_grf
    ) 
  )
  
  tibble::rownames_to_column(res, "estimator")
  
  return(res)
}                                                                               

### RUNNING THE FUNCTION
system.time(
  results <- lapply(
    1:number_sims, 
    function(x) cluster_sim(nsim=x,sample_size=sample_size)
    )
  )
results <- do.call(rbind, results)

write_csv(results, file = here("data","simulation_results-phase1.csv"))