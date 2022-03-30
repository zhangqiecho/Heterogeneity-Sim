#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## NOTE: We can run this program locally from the command line, and thus leave the package 
## installation as is. Once we run on cluster, need to change library location.

## where is the library?
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","skimr","here","mvtnorm","latex2exp","earth",
   "readxl","VGAM", "ranger","xgboost","mgcv","glmnet","NbClust","factoextra",
  "SuperLearner", "AIPW", "dplyr", "cluster", "ggplot2")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

## FUNCTION SET-UP
## TRUE VALUES
true1 <- 6
true2 <- 3

number_sims <- as.numeric(args[1])

cat(paste("Number of Simulations:", number_sims, "\n"))

sample_size <- as.numeric(args[2])

cat(paste("Sample Size for Each Sim:", sample_size), "\n")

random_argument_for_checking <- args[3]

# this should say "hello_world!"
cat(paste("Do you have anything to say?:",random_argument_for_checking), "\n")

## FUNCTION
cluster_sim <- function(nsim, sample_size){
  
  n = sample_size
  p = 5
  set.seed(nsim)
  
  ## CONFOUNDERS (C = 5)
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  c     <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESIGN MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  parms3 <- rep(3,5)
  parms4 <- rep(log(1.5),5)
  beta   <- parms3; beta <- c(120, beta)
  
  # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
  piMatT  <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  piMatT2 <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  theta   <- c(-.5,parms4)
  theta2  <- c(-.5,parms4)
  mu      <- muMatT%*%beta
  
  # PROPENSITY SCORE MODEL
  ## pi is an actual parameter in base R, so good practice not to overwrite
  pi_x   <- expit(piMatT%*%theta)
  pi_m <- expit(piMatT2%*%theta2)
  x    <- rbinom(n,1,pi_x)
  m    <- rbinom(n,1,pi_m)
  
  # OUTCOME MODEL: EXPOSURE VALUE UNDER M == 0 IS 6; VALUE UNDER M == 1 IS 3
  y0 <- 86 + 3*m[x==0] + 4*c[x==0,1] + 3*c[x==0,2] + 2*c[x==0,3] + 1*c[x==0,4] + 1*c[x==0,5] + rnorm(sum(1-x),0,20)                        #the difference in these two models may be explaining the quasi/complete separation
  y1 <- 88 + 6*m[x==1] + 4*c[x==1,1] + 3*c[x==1,2] + 2*c[x==1,3] + 1*c[x==1,4] + 1*c[x==1,5] + rnorm(sum(x),0,20)
  
  #y0 <- 88 + 3*m[x==0] + 4*c[x==0,1] + 3*c[x==0,2] + 2*c[x==0,3] + 1*c[x==0,4] + 1*c[x==0,5] + rnorm(sum(1-x),0,20)                        #the difference in these two models may be explaining the quasi/complete separation
  #y1 <- 86 + 6*m[x==1] + 1*c[x==1,1] + 1*c[x==1,2] - 1.5*c[x==1,3] - 5*c[x==1,4] + 1*c[x==1,5] + rnorm(sum(x),0,15)
  #ATE: y1-y0 --> average over all of these differences in effects (rather than pure effect measure modification)
  d0 <- data.frame(y0,m[x==0],c[x==0,],x=0);names(d0) <- c("y","m",paste0("c",1:5),"x")
  d1 <- data.frame(y1,m[x==1],c[x==1,],x=1);names(d1) <- c("y","m",paste0("c",1:5),"x")
  
  # DATA
  dat <- tibble(rbind(d0,d1))
  dat
  
  ## DATA GENERATION FINISHED
  
  ## DATA ANALYSIS START
  
  ## EXAMPLE 1
  ## simplest approach (for illustration), not recommended in practice
  mod1_x1 <- glm(y ~ m + c1 + c2 + c3 + c4 + c5, data=subset(dat,x==1), family=gaussian(link="identity"))
  mod1_x0 <- glm(y ~ m + c1 + c2 + c3 + c4 + c5, data=subset(dat,x==0), family=gaussian(link="identity"))
  summary(mod1_x1)
  summary(mod1_x0)
  
  # dat1 <- transform(dat,x=1)
  # dat0 <- transform(dat,x=0)
  # 
  # head(dat1)
  # head(dat0)
  
  dat$mu1 <- predict(mod1_x1,dat,type="response")
  dat$mu0 <- predict(mod1_x0,dat,type="response")
  
  mean(dat$mu1-dat$mu0)
  
  clust_dat <- dat %>% select(mu0, mu1)
  
  clust_dat
  
  # now search for clusters
  #trying the NbClust package
  pdf(file = NULL)
  clust <- NbClust(data = clust_dat, diss=NULL, distance = "euclidean", min.nc = 2, max.nc = 15, method = "kmeans", index = "all", alphaBeale = 0.1)
  dev.off()

  df <- data.frame(clust$Best.nc[1,1:26])
  names(df) <- "number_clust"
  row.names(df) <- NULL
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  cluster_num <- getmode(df$number_clust)
  kmeans_res <- kmeans(clust_dat, centers = cluster_num, max.iter=100, nstart = 25)
  
  dat$cluster <- kmeans_res$cluster # check: this should give us a cluster number for each obs in the original data
  dat$diff <- dat$mu1-dat$mu0
  cluster_effects <- aggregate(dat$diff, list(dat$cluster), mean) #we should technically see that it is yielding two clusters and the average output to be equal to the truth (if working appropriately)
  
  if(cluster_num > 2){
    cluster_reg <- vglm(cluster~.,data=subset(dat,select=-c(y,mu0,mu1,diff)),family="multinomial") 
    res4 <- c(summary(cluster_reg)@coefficients, diag(sqrt(summary(cluster_reg)@cov.unscaled)))
    # check to see that subset gives a dataset with only x, m, c, and cluster
  } else{
    cluster_reg <- glm(as.numeric(cluster==2)~., data=subset(dat, select=-c(y,mu0,mu1,diff)),family=binomial(link="logit"))
    res4 <- c(summary(cluster_reg)$coefficients[,1:2])
  }
  
  res1 <- c(mod1_x0$coefficients)
  res2 <- c(mod1_x1$coefficients)
  res3 <- c(cluster_effects)


# check: we want point estimates and SE
  res <- list(res1,res2,res3,res4)
  return(res)
}                                                                               

### RUNNING THE FUNCTION
system.time(results <- lapply(1:number_sims, function(x) cluster_sim(nsim=x,sample_size=sample_size)))

saveRDS(results, file = here("data","Results.Rds"))

print(results[[1]])