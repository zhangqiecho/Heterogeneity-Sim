#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## NOTE: We can run this program locally from the command line, and thus leave the package 
## installation as is. Once we run on cluster, need to change library location.

start.time <- Sys.time()

## where is the library?
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","skimr",
              "here","mvtnorm","latex2exp","earth",
              "readxl", "ranger","xgboost",
              "mgcv","glmnet","NbClust","factoextra",
              "dplyr", "cluster", "foreach","doParallel",
              "ggplot2", "estimatr","grf", "DoubleML",
              "mlr3","mlr3learners","mlr3pipelines")

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
remotes::install_github('tlverse/tlverse')
library(tmle3)
library(sl3)
remotes::install_github('yqzhong7/AIPW')
library(AIPW)

#args = c(2,1000)

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

number_sims <- as.numeric(args[1])
cat(paste("Number of Simulations:", number_sims, "\n"))
sample_size <- as.numeric(args[2])
cat(paste("Sample Size for Each Sim:", sample_size), "\n")

## FUNCTION
cluster_sim <- function(nsim, sample_size, param1, param2){

  param1 <- exp(param1)
  param2 <- exp(param2)
  
  n = sample_size
  p = 5
  set.seed(nsim)
  print(paste("simulation number",nsim))
  
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
  mu_y <- expit(-log((1/.2)-1) - log(param1)*pi_x - log(7.5)*pi_m + 
                  log(param2)*pi_x*pi_m + log(param1)*x+ log(7.5)*m - log(param2)*x*m)
  
  y <- rbinom(n, 1, mu_y)

  dat <- data.frame(y,m,c,x)
  names(dat) <- c("y","m",paste0("c",1:5),"x")
  
  ## DATA GENERATION FINISHED; DATA ANALYSIS START
  
  # TRUE MODEL
  mod_true <- lm(y ~ x + m + x*m + c1 + c2 + c3 + c4 + c5, data = dat)
  c_matrix1 <- c(0,1,0,0,0,0,0,0,0) # ATE among m = 0
  c_matrix2 <- c(0,1,0,0,0,0,0,0,1) # ATE among m = 1; NB: r puts interaction coefficient at end
  
  ate_m0 <- c(coef(mod_true) %*% c_matrix1, 
              sqrt(t(c_matrix1) %*% vcov(mod_true) %*% c_matrix1))
  
  ate_m1 <- c(coef(mod_true) %*% c_matrix2, 
              sqrt(t(c_matrix2) %*% vcov(mod_true) %*% c_matrix2))
  
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
  
  ate_m1_grf <- average_treatment_effect(cf, subset = m == 1) 
  ate_m0_grf <- average_treatment_effect(cf, subset = m == 0) 
  
  # TMLE & AIPW
  
  print(paste("Constructing learners for TMLE & AIPW ... "))
  
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_glm <- make_learner(Lrnr_glm)
  
  # ranger learner
  grid_params <- list(num.trees = c(250, 1000),
                      mtry = c(3,4),
                      min.node.size = c(50,100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  lrnr_ranger <- vector("list", length = nrow(grid))
  for(i in 1:nrow(grid)){
    lrnr_ranger[[i]] <- make_learner(Lrnr_ranger, 
                                     num.trees=grid[i,]$num.trees, 
                                     mtry=grid[i,]$mtry,
                                     min.node.size=grid[i,]$min.node.size)
  }
  ########################################################################
  #lrnr_ranger <- make_learner(Lrnr_ranger)  ###########  FLAG! ###########
  ########################################################################
  
  # glmnet learner
  grid_params <- seq(0,1,by=.25)
  lrnr_glmnet <- vector("list", length = length(grid_params))
  for(i in 1:length(grid_params)){
    lrnr_glmnet[[i]] <- make_learner(Lrnr_glmnet, alpha = grid_params[i])
  }
  ########################################################################
  #lrnr_glmnet <- make_learner(Lrnr_glmnet)  ###########  FLAG! ###########
  ########################################################################
  
  # xgboost learner
  grid_params <- list(max_depth = c(4, 6),
                      eta = c(0.1, 0.2),
                      nrounds = c(100, 500)
  )
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  lrnr_xgboost <- vector("list", length = nrow(grid))
  for(i in 1:nrow(grid)){
    lrnr_xgboost[[i]] <- make_learner(Lrnr_xgboost, max_depth=grid[i,]$max_depth, eta=grid[i,]$eta)
  }
  ########################################################################
  #lrnr_xgboost <- make_learner(Lrnr_xgboost)  ##########  FLAG! ##########
  ########################################################################
  
  # earth learner
  grid_params <- c(3,4,5)
  lrnr_earth <- vector("list", length = length(grid_params))
  for(i in 1:length(grid_params)){
    lrnr_earth[[i]] <- make_learner(Lrnr_earth, degree = grid_params[i])
  }
  ########################################################################
  ##lrnr_earth <- make_learner(Lrnr_earth)  ############  FLAG! ############
  ########################################################################
  
  sl_ranger <- make_learner(Stack, unlist(list(lrnr_ranger), 
                                    recursive = TRUE))
  sl_xgboost <- make_learner(Stack, unlist(list(lrnr_xgboost), 
                                          recursive = TRUE))
  sl_tree <- make_learner(Stack, unlist(list(lrnr_ranger,
                                               lrnr_xgboost), 
                                          recursive = TRUE))
  sl_ <- make_learner(Stack, unlist(list(lrnr_ranger,
                                         lrnr_xgboost,
                                         lrnr_mean,
                                         lrnr_glm,
                                         lrnr_earth,
                                         lrnr_glmnet), 
                                        recursive = TRUE))
  
  # DEFINE SL_Y AND SL_A 
  # We only need one, because they're the same
  
  learners_ <- c(sl_ranger, sl_xgboost, sl_tree, sl_)
  tmle0_list <- tmle1_list <- aipw0_list <- aipw1_list <- list()

  for (i in 1:length(learners_)){
    
  Q_learner <- Lrnr_sl$new(learners = learners_[i], 
                           metalearner = Lrnr_nnls$new(convex=T))
  g_learner <- Lrnr_sl$new(learners = learners_[i], 
                           metalearner = Lrnr_nnls$new(convex=T))
  learner_list <- list(Y = Q_learner,
                       A = g_learner)
  
  ######################################################################
  
  # PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
  ate_spec <- tmle_OR(contrast_level = 1, baseline_level = 0)
  
  nodes_ <- list(W = c("c1", "c2", "c3", "c4", "c5"), 
                 A = "x", 
                 Y = "y")
  
  # RUN TMLE3 
  print(paste("TMLE Version",i))
  set.seed(123)
  tmle0_list[[i]] <- tmle3(ate_spec, 
                           subset(dat, m==0, select = -m), 
                           nodes_, 
                           learner_list)
  tmle1_list[[i]] <- tmle3(ate_spec, 
                           subset(dat, m==1, select = -m), 
                           nodes_, 
                           learner_list)

  print(paste("AIPW Version",i))
  aipw0_list[[i]] <- AIPW_tmle$new(A=dat[dat$m==0,]$x,
                                   Y=dat[dat$m==0,]$y,
                                   tmle_fit = tmle0_list[[i]],
                                   verbose = T)$summary()
  
  aipw1_list[[i]] <- AIPW_tmle$new(A=dat[dat$m==1,]$x,
                                   Y=dat[dat$m==1,]$y,
                                   tmle_fit = tmle1_list[[i]],
                                   verbose = TRUE)$summary()
  
  }
  
  ate_m0_tmle_ranger <- cbind(tmle0_list[[1]]$estimates[[2]]$psi - tmle0_list[[1]]$estimates[[1]]$psi, 
                              sqrt(var(tmle0_list[[1]]$estimates[[2]]$IC - tmle0_list[[1]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m1_tmle_ranger <- cbind(tmle1_list[[1]]$estimates[[2]]$psi - tmle1_list[[1]]$estimates[[1]]$psi, 
                              sqrt(var(tmle1_list[[1]]$estimates[[2]]$IC - tmle1_list[[1]]$estimates[[1]]$IC)/sum(1-m)))

  ate_m0_tmle_xgboost <- cbind(tmle0_list[[2]]$estimates[[2]]$psi - tmle0_list[[2]]$estimates[[1]]$psi, 
                               sqrt(var(tmle0_list[[2]]$estimates[[2]]$IC - tmle0_list[[2]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m1_tmle_xgboost <- cbind(tmle1_list[[2]]$estimates[[2]]$psi - tmle1_list[[2]]$estimates[[1]]$psi, 
                               sqrt(var(tmle1_list[[2]]$estimates[[2]]$IC - tmle1_list[[2]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m0_tmle_tree <- cbind(tmle0_list[[3]]$estimates[[2]]$psi - tmle0_list[[3]]$estimates[[1]]$psi, 
                            sqrt(var(tmle0_list[[3]]$estimates[[2]]$IC - tmle0_list[[3]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m1_tmle_tree <- cbind(tmle1_list[[3]]$estimates[[2]]$psi - tmle1_list[[3]]$estimates[[1]]$psi, 
                            sqrt(var(tmle1_list[[3]]$estimates[[2]]$IC - tmle1_list[[3]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m0_tmle_stacked <- cbind(tmle0_list[[4]]$estimates[[2]]$psi - tmle0_list[[4]]$estimates[[1]]$psi, 
                               sqrt(var(tmle0_list[[4]]$estimates[[2]]$IC - tmle0_list[[4]]$estimates[[1]]$IC)/sum(1-m)))
  
  ate_m1_tmle_stacked <- cbind(tmle1_list[[4]]$estimates[[2]]$psi - tmle1_list[[4]]$estimates[[1]]$psi, 
                               sqrt(var(tmle1_list[[4]]$estimates[[2]]$IC - tmle1_list[[4]]$estimates[[1]]$IC)/sum(1-m)))
  
  
  ate_m0_aipw_ranger <- aipw0_list[[1]]$result["Risk Difference",1:2]
  
  ate_m1_aipw_ranger <- aipw1_list[[1]]$result["Risk Difference",1:2]
  
  ate_m0_aipw_xgboost <- aipw0_list[[2]]$result["Risk Difference",1:2]
  
  ate_m1_aipw_xgboost <- aipw1_list[[2]]$result["Risk Difference",1:2]
  
  ate_m0_aipw_tree <- aipw0_list[[3]]$result["Risk Difference",1:2]
    
  ate_m1_aipw_tree <- aipw0_list[[3]]$result["Risk Difference",1:2]
    
  ate_m0_aipw_stacked <- aipw0_list[[4]]$result["Risk Difference",1:2]
  
  ate_m1_aipw_stacked <- aipw0_list[[4]]$result["Risk Difference",1:2]
  
  print(paste("Creating res object ... "))
  
  res <- data.frame(
	estimator = c("CATE_m0_glm","CATE_m1_glm",
	              "CATE_m0_grf","CATE_m1_grf",
	              "CATE_m0_tmle_ranger","CATE_m1_tmle_ranger",
	              "CATE_m0_tmle_xgboost","CATE_m1_tmle_xgboost",
	              "CATE_m0_tmle_tree","CATE_m1_tmle_tree",
	              "CATE_m0_tmle_stacked","CATE_m1_tmle_stacked",
	              "CATE_m0_aipw_ranger", "CATE_m1_aipw_ranger",
	              "CATE_m0_aipw_xgboost","CATE_m1_aipw_xgboost",
	              "CATE_m0_aipw_tree",   "CATE_m1_aipw_tree",
	              "CATE_m0_aipw_stacked","CATE_m1_aipw_stacked"),
    rbind(
      ate_m0,ate_m1,
      ate_m0_grf,ate_m1_grf,
      ate_m0_tmle_ranger,ate_m1_tmle_ranger,
      ate_m0_tmle_xgboost,ate_m1_tmle_xgboost,
      ate_m0_tmle_tree,ate_m1_tmle_tree,
      ate_m0_tmle_stacked,ate_m1_tmle_stacked,
      ate_m0_aipw_ranger,ate_m1_aipw_ranger,
      ate_m0_aipw_xgboost,ate_m1_aipw_xgboost,
      ate_m0_aipw_tree,ate_m1_aipw_tree,
      ate_m0_aipw_stacked,ate_m1_aipw_stacked
    ),
	parameter1 = param1, 
	parameter2 = param2,
	simulation_number = nsim
  )

  print(paste("Writing to file ... "))
  
  write_csv(res, file = here("data",paste("simulation_results_",nsim,".csv")), append = T)
  
  print(paste("simulation", nsim, "end"))
  return(res)
}                                  

### RUNNING THE FUNCTION
phi_dat <- list(phi = seq(-4,4,by=2), 
                phi2 = seq(-6,6,by=2))

data_grid <- expand.grid(phi_dat, KEEP.OUT.ATTRS = FALSE)

parallel::detectCores()

n.cores <- parallel::detectCores()

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()

par_res <- foreach(i = 1:nrow(data_grid)) %dopar% {
#par_res <- foreach(i = 1:n.cores) %dopar% {
  
  print(paste("for loop iteration",i))
  
  results <- do.call(rbind,
                      lapply(1:number_sims, 
                                function(x) cluster_sim(nsim=x,
                                                        sample_size=sample_size, 
                                                        param1 = data_grid[i,]$phi, 
                                                        param2 = data_grid[i,]$phi2)
                             )
                      )
  
}

parallel::stopCluster(cl = my.cluster)

par_res <- do.call(rbind, par_res)

rownames(par_res) <- NULL

write_rds(par_res, 
          file = here("data",paste("simulation_results_",number_sims,".rds"))
          )

end.time <- Sys.time()

duration_time <- end.time - start.time

dur_dat <- data.frame(run_time = duration_time,
                      number_sims = number_sims)

write_csv(dur_dat, file = here("data", paste0("run_time_",number_sims,".csv")))