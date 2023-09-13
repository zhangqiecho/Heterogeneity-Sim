#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## NOTE: We can run this program locally from the command line, and thus leave the package 
## installation as is. Once we run on cluster, need to change library location.

start.time <- Sys.time()

## where is the library?
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","SuperLearner","boot",
              "here","mvtnorm","foreach","doParallel","broom",
              "ggplot2","grf", "lmtest", "sandwich")

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

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

number_sims <- as.numeric(args[1])
cat(paste("Number of Simulations:", number_sims, "\n"))
sample_size <- as.numeric(args[2])
cat(paste("Sample Size for Each Sim:", sample_size), "\n")

## FUNCTION
cluster_sim <- function(nsim, sample_size, param1, param2, c_number){

  set.seed(nsim)
  ## FOCUS ON PROBLEM OF USING CF AND DR LEARNER TO FIND THE RIGHT MODIFIER
  ## RE-PROGRAM TO THAT P CAN BE A VARIABLE IN THE SIMS
  ## REMOVE COMPLEXITY IN SUPER LEARNER, JUST USE THIS: sl_
  
  param1 <- exp(param1) # these are huge effects. need to change...
  param2 <- exp(param2)
  
  n = sample_size
  p = c_number
  print(paste("simulation number",nsim))
  
  ## CONFOUNDERS
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  c     <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESIGN MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  parmsC <- rep(log(1.5),c_number)
  
  # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
  piMatX  <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  piMatM  <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  thetaX  <- c(-.5,parmsC) # assumed coefficients of propensity score model
  thetaM  <- c(-.5,parmsC)
  
  # PROPENSITY SCORE MODEL
  ## pi is an actual parameter in base R, so good practice not to overwrite
  pi_x <- expit(piMatX%*%thetaX)
  pi_m <- expit(piMatM%*%thetaM)
  x    <- rbinom(n,1,pi_x)
  m    <- rbinom(n,1,pi_m)
  
  # OUTCOME MODEL: EXPOSURE VALUE UNDER M == 0 IS 7.5; VALUE UNDER M = 1 IS -7.5
  mu_y <- expit(-log((1/.2)-1) - log(param1)*pi_x - log(7.5)*pi_m + log(param2)*pi_x*pi_m + # this is an approximation to the balanced intercept (c's are not included)
                 log(param1)*x+ log(7.5)*m - log(param2)*x*m + # this is the exposure and modifer effect
                  c%*%parmsC # this is the confounder-outcome part
                )
  
  y <- rbinom(n, 1, mu_y)

  dat <- data.frame(y,m,c,x)
  names(dat) <- c("y","m",paste0("c",1:c_number),"x")
  
  ## DATA GENERATION FINISHED; DATA ANALYSIS START

  # ORACLE MODEL, G COMP
  mod_form <- as.formula(paste("y ~ x + m + x*m +", paste0("c", 1:c_number, collapse = "+")))
  mod_true <- glm(mod_form, data = dat, family = binomial(link = "logit"))
  
  # c_matrix1 <- c(0, (as.matrix(attr(summary(mod_true)$terms,"term.labels")) == c("x")))
  # c_matrix2 <- c_matrix1 + c(0, (as.matrix(attr(summary(mod_true)$terms,"term.labels")) == c("x:m")))
  # 
  # ate_m0 <- c(coef(mod_true) %*% c_matrix1, 
  #             sqrt(t(c_matrix1) %*% vcov(mod_true) %*% c_matrix1))
  # 
  # ate_m1 <- c(coef(mod_true) %*% c_matrix2, 
  #             sqrt(t(c_matrix2) %*% vcov(mod_true) %*% c_matrix2))
  
  mu1 <- mean(predict(mod_true, newdata = transform(dat, x = 1), type = "response"))
  mu0 <- mean(predict(mod_true, newdata = transform(dat, x = 0), type = "response"))
  
  ate_gcomp <- mu1 - mu0
  
  mu11 <- mean(predict(mod_true, newdata = transform(dat, x = 1, m = 1), type = "response"))
  mu01 <- mean(predict(mod_true, newdata = transform(dat, x = 0, m = 1), type = "response"))
  
  ate_m1_gcomp <- mu11 - mu01
  
  mu10 <- mean(predict(mod_true, newdata = transform(dat, x = 1, m = 0), type = "response"))
  mu00 <- mean(predict(mod_true, newdata = transform(dat, x = 0, m = 0), type = "response"))
  
  ate_m0_gcomp <- mu10 - mu00
  
  boot_func <- function(data, index){
    
    boot_dat <- data[index, ]
    
    mod_true_ <- glm(mod_form, data = boot_dat, family = binomial(link = "logit"))
    
    mu1_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1), type = "response"))
    mu0_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0), type = "response"))
    
    ate_gcomp_ <- mu1_ - mu0_
    
    mu11_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1, m = 1), type = "response"))
    mu01_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0, m = 1), type = "response"))
    
    ate_m1_gcomp_ <- mu11_ - mu01_
    
    mu10_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1, m = 0), type = "response"))
    mu00_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0, m = 0), type = "response"))
    
    ate_m0_gcomp_ <- mu10_ - mu00_

    return(c(ate_gcomp_,
             ate_m0_gcomp_,
             ate_m1_gcomp_))
        
  }
  
  boot_res <- boot(boot_func, data = dat, R = 500)
  
  ate_gcomp <- c(ate_gcomp, sd(boot_res$t[,1]))
  names(ate_gcomp) <- c("estimate", "std.err")
  
  ate_m0_gcomp <- c(ate_m0_gcomp, sd(boot_res$t[,2]))
  names(ate_m0_gcomp) <- c("estimate", "std.err")
  
  ate_m1_gcomp <- c(ate_m1_gcomp, sd(boot_res$t[,3]))
  names(ate_m1_gcomp) <- c("estimate", "std.err")

  ## HERE WE HAVE CAUSAL FOREST CODE TO IDENTIFY APPROPORATE MODIFERS
  ## 
  covariates_matrix <- dat[,c("m", paste0("c", 1:c_number))]
  covariates_matrix_x <- dat[,c("x", "m", paste0("c", 1:c_number))]
  outcome <-  as.matrix(dat$y)
  exposure <- as.matrix(dat$x)

  # construct a sample size variable 
  n <- nrow(dat)
  
  # We use 10-fold cross fitting, and determine precisely which of the observations get assigned to each fold:
  num.folds <- 10
  folds <- sort(seq(n) %% num.folds) + 1
  
  # we construct a causal forest using the `causal_forest` function from the grf package:
  forest <- causal_forest(X = covariates_matrix, 
                          Y = outcome, 
                          W = exposure, 
                          num.trees = 2000,
                          honesty = TRUE,
                          min.node.size = 10,
                          alpha = .05, 
                          imbalance.penalty = 0,
                          stabilize.splits = TRUE,
                          tune.parameters = "all",
                          tune.num.trees = 200,
                          tune.num.reps = 50,
                          tune.num.draws = 1000,
                          compute.oob.predictions = TRUE,
                          num.threads = 10,
                          seed = 123,
                          clusters = folds)
  
  tau.hat = predict(forest)$predictions
  
  ## COMPUTE THE ATE FROM THE CAUSAL FOREST ALGORITHM
  # we then use the `average_treatment_effect()` function to obtain an ATE estimate
  cf_ate <- average_treatment_effect(forest)
  
  blp_cf <- tidy(best_linear_projection(forest, covariates_matrix)) %>% 
    filter(term != "(Intercept)") %>% 
    arrange(desc(abs(statistic))) %>% 
    mutate(parameter1 = param1,
           parameter2 = param2, 
           sample_size = n, 
           confounder_number = c_number, 
           seed_number = nsim, 
           method = "BLP_causal_forest")
  ## need to pick the "strong" ones from this, and store the results
  ## what we should do is rank the |test statistic| and pick the first one to check if it's m
  ## maybe keep this whole thing??
  
  aipw_scores <- get_scores(forest,
                           subset = NULL,
                           debiasing.weights = NULL,
                           num.trees.for.weights = 2000
                           )
  
  fmla <- as.formula(aipw_scores ~ factor(m))
  score_reg <- lm(fmla, 
                  data=transform(dat, aipw_scores=aipw_scores))
  
  c_matrix1 <- c(1, 0)
  c_matrix2 <- c(1, 1)

  cf_cate_m0 <- data.frame(t(c(coef(score_reg) %*% c_matrix1, 
                             sqrt(t(c_matrix1) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix1))
                             ))
  names(cf_cate_m0) <- c("estimate", "std.err")

  cf_cate_m1 <- data.frame(t(c(coef(score_reg) %*% c_matrix2, 
                               sqrt(t(c_matrix2) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix2))
  ))
  names(cf_cate_m1) <- c("estimate", "std.err")
  

  
  # DO THE EXACT SAME FOR DRLEARNER.
  
  
  mean_learner <- "SL.mean"
  glm_learner <- "SL.glm"
  
  # ranger learner
  ranger_learner <- "SL.ranger"
  
  # glmnet learner
  glmnet_learner <- "SL.glmnet"
  
  # xgboost learner
  xgboost_learner <- "SL.xgboost"
  
  # earth learner
  earth_learner <- "SL.earth"
  
  sl_lib <- c(ranger_learner,
              glmnet_learner, 
              xgboost_learner, 
              mean_learner,
              glm_learner)
  
  # Specify the number of folds for V-fold cross-validation
  # Use same folds as used for causal_forest function
  # Doing cross-validation this way automatically deploys cross-fitting
  fold_dat <- tibble(id = 1:n, folds)
  fold_index <- split(fold_dat$id, fold_dat$folds)
  
  fit_mu <- CV.SuperLearner(Y = outcome,
                            X = covariates_matrix_x, 
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),
                            control = list(saveCVFitLibrary = T),
                            parallel = "seq",
                            verbose = T)
  
  fit_pi <- CV.SuperLearner(Y = exposure,
                            X = covariates_matrix,
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),#, stratifyCV = TRUE),
                            control = list(saveCVFitLibrary = T),
                            parallel = "seq",
                            verbose = T)
  
  pscore <- as.matrix(fit_pi$SL.predict)
  
  mu_hat <- as.matrix(fit_mu$SL.predict)
  
  mu_hat1 <- NULL
  for(i in 1:num.folds){
    mu_hat1 <- rbind(mu_hat1, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_matrix_x[fold_index[[i]],], x = 1), 
                             onlySL=T)$pred)
  }
  
  mu_hat0 <- NULL
  for(i in 1:num.folds){
    mu_hat0 <- rbind(mu_hat0, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_matrix_x[fold_index[[i]],], x = 0), 
                             onlySL=T)$pred)
  }
  
  ## aipw
  aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1){
    aipw_score <- ((2*exposure - 1)*(outcome - mu_hat))/((2*exposure - 1)*pscore + (1 - exposure)) + (mu_hat1 - mu_hat0)
    return(aipw_score)
  }
  
  aipw_score <- aipw_func(exposure, 
                          outcome, 
                          pscore, 
                          mu_hat, 
                          mu_hat0, 
                          mu_hat1)
  colnames(aipw_score) <- NULL
  
  aipw_psi <- mean(aipw_score)
  
  aipw_se <- sd(aipw_score)/sqrt(n)
  
  aipw_ate <- data.frame(t(c(aipw_psi, aipw_se)))
  names(aipw_ate) <- c("estimate", "std.err")
  
  ## best linear projection
  score_dat <- data.frame(aipw_score, covariates_matrix)
  score_fit <- lm(aipw_score ~ ., data = score_dat)
  
  blp_aipw <- coeftest(score_fit, vcov = vcovHC(score_fit, type = "HC3"))
  
  blp_drl <- tidy(blp_aipw) %>% 
    filter(term != "(Intercept)") %>% 
    arrange(desc(abs(statistic))) %>% 
    mutate(parameter1 = param1,
           parameter2 = param2, 
           sample_size = n, 
           confounder_number = c_number,
           seed_number = nsim, 
           method = "BLP_dr_learner")
  
  fmla <- as.formula(aipw_score ~ factor(m))
  score_reg <- lm(fmla, 
                  data=transform(dat, aipw_scores=aipw_score))
  
  c_matrix1 <- c(1, 0)
  c_matrix2 <- c(1, 1)
  
  drl_cate_m0 <- data.frame(t(c(coef(score_reg) %*% c_matrix1, 
                               sqrt(t(c_matrix1) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix1))
  ))
  names(drl_cate_m0) <- c("estimate", "std.err")
  
  drl_cate_m1 <- data.frame(t(c(coef(score_reg) %*% c_matrix2, 
                               sqrt(t(c_matrix2) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix2))
  ))
  names(drl_cate_m1) <- c("estimate", "std.err")

  # objects to save:
  
  res <- rbind(
    data.frame(estimator = "gComp_ate",            seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, t(ate_gcomp)),
    data.frame(estimator = "gComp_cate_m0",        seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, t(ate_m0_gcomp)),
    data.frame(estimator = "gComp_cate_m1",        seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, t(ate_m1_gcomp)),
    data.frame(estimator = "causal_forest_ate",    seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, t(cf_ate)),
    data.frame(estimator = "causal_forest_ate_m0", seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, cf_cate_m0),
    data.frame(estimator = "causal_forest_ate_m1", seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, cf_cate_m1),
    data.frame(estimator = "DRLearner_ate",        seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, aipw_ate),
    data.frame(estimator = "DRLearner_ate_m0",     seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, drl_cate_m0),
    data.frame(estimator = "DRLearner_ate_m1",     seed_nubmer = nsim, sample_size = n, parameter1 = param1, parameter2 = param2, drl_cate_m1)
  )
  
  func_output <- list(blp_cf,
                      blp_drl,
                      res)
  
  return(func_output)

}

cnum=5
number_sims=3
sample_size=1000
run_loop=1

results <- lapply(1:number_sims, 
                  function(x) cluster_sim(nsim = x,
                                          c_number = cnum,
                                          sample_size = sample_size, 
                                          param1 = 2, 
                                          param2 = 2)
)

### RUNNING THE FUNCTION
phi_dat <- list(phi = seq(-2,2,by=2), 
                phi2 = seq(-2,2,by=2))

data_grid <- expand.grid(phi_dat, KEEP.OUT.ATTRS = FALSE)

parallel::detectCores()

n.cores <- parallel::detectCores() - 2

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()

number_sims <- as.numeric(args[1])
cat(paste("Number of Simulations:", number_sims, "\n"))
sample_size <- as.numeric(args[2])
cat(paste("Sample Size for Each Sim:", sample_size), "\n")
run_loop <- as.numeric(args[3])
cat(paste("The Loop Currently Running is:", run_loop), "\n")
cnum <- as.numeric(args[4])
cat(paste("The Number of Confounders is:", cnum), "\n")

cnum=5
number_sims=3
sample_size=5000
run_loop=1


par_res <- foreach(i = 1:nrow(data_grid)) %dopar% {
  
  print(paste("for loop iteration",i))
  
  results <- lapply(1:number_sims, 
                    function(x) cluster_sim(nsim = x,
                                            c_number = cnum,
                                            sample_size = sample_size, 
                                            param1 = data_grid[i,]$phi, 
                                            param2 = data_grid[i,]$phi2)
                             )
  
}

parallel::stopCluster(cl = my.cluster)

## first bracket: indexing over parameters phi and phi2
## second bracket: indexing over simulation numbers
## third bracket: indexing over result - blp_cf, blp_drl, and simulation results

k <- 1:(nrow(data_grid)*number_sims)
blp_cf <- blp_drl <- sim_res <- NULL
for(i in 1:nrow(data_grid)){
  for(j in 1:number_sims){
    blp_cf[k]  <- list(par_res[[i]][[j]][[1]])
    blp_drl[k] <- list(blp_drl, par_res[[i]][[j]][[2]])
    sim_res[k] <- list(sim_res, par_res[[i]][[j]][[3]])
  }
}

blp_cf <- do.call(rbind,blp_cf)

blp_drl <- par_res[[1]][[1]][[2]]

sim_res <- par_res[[2]][[2]][[3]]

write_rds(par_res, 
          file = here("data",paste0("simulation_results_",number_sims,".rds"))
          )

end.time <- Sys.time()

duration_time <- end.time - start.time

dur_dat <- data.frame(run_time = duration_time,
                      number_sims = number_sims)

write_csv(dur_dat, 
          file = here("data", 
                      paste0("run_time_",number_sims,run_loop,".csv")))