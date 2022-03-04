###################################################################################################################
#
# Purpose: Simulation Model 
#
# Model 1: Y <- intercept + x*a + m*b - d*x*m + c1...c5 + error
#
# Last Update: 29 Jan 2022
#
###################################################################################################################

####################################################
#               LOAD IN PACKAGES
####################################################

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

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)


## FUNCTION SET-UP
## TRUE VALUES
true1 <- 6
true2 <- 3

#nsim <- 5                                                                        #cancelled out when function is ran

#DUMMY TABLE TO EXPORT THE RESULTS OF THE SIMULATION
#cols <- c("point_est1", "point_est2", "num_clust")
#cols <- c(cols, "nsim")
#res.est <- data.frame(matrix(nrow=nsim,ncol=length(cols)))
#colnames(res.est) <- cols
#res.se <- res.est

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

## FUNCTION
cluster_sim <- function(nsim, sample_size){                             
  i = nsim
  #i = 6
  n = sample_size                                                                      
  #n = 100
  p = 5
  set.seed(i)                                                                 
  
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
  pi   <- expit(piMatT%*%theta)
  pi_m <- expit(piMatT2%*%theta2)
  x    <- rbinom(n,1,pi)
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
  clust <- NbClust(data = clust_dat, diss=NULL, distance = "euclidean", min.nc = 2, max.nc = 15, method = "kmeans", index = "all", alphaBeale = 0.1)
  
  df <- data.frame(clust$Best.nc[1,1:26])
  names(df) <- "number_clust"
  row.names(df) <- NULL
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  cluster_num <- getmode(df$number_clust)
  kmeans_res <- kmeans(clust_dat, centers = cluster_num, nstart = 25)
  
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

summary(cluster_reg)
#cluster_sim(nsim = 6, sample = 100)
##############################################################

### RUNNING THE FUNCTION
system.time(results <- lapply(1:1000, function(x) cluster_sim(nsim=x,sample_size=1000))) #Error: $ operator not defined for this S4 class [instead of $, use]
saveRDS(results, file = 'Results.Rds')

results[[1]] #output from the first simulation run
#estimated effect of m is [[2]]-[[1]] for the 
####
#
# Next Steps: 
#
# 10 runs with N = 1,000 (338.69 seconds) 5 mins
# 20 runs with N = 1,000 (667.31 seconds) 10 mins
# 40 runs with N = 1,000 (1300.45 seconds) 20 mins
# 60 runs with N = 1,000 (2420.32 seconds) 40 mins
# 
#
runs <- c(10, 20, 40, 60)
minutes <- c(5, 10, 20, 40)
sim_time <- data.frame(runs, minutes)
ggplot(sim_time, aes(runs, minutes)) +
  geom_point() +
  ggtitle("Time to run simulations with sample size of 1,000")
#
#
####
#

### COMPILING THE RESULTS
res <- do.call(rbind, results)
res <- as.data.frame(res)

#What percentage of runs result in n = 2 clusters?
#select the number of clusters from each run
number_of_cluster <-matrix(NA, ncol = 1)
clusters <- for(i in 1:1000){
  #number_of_clusters <- rbind(number_of_clusters, max(results[[i]][[3]]$Group.1))
  number_of_cluster <- rbind(number_of_cluster,as.matrix(max(results[[i]][[3]]$Group.1)))
}
number_of_cluster <- number_of_cluster[-1,]
table(number_of_cluster)

saveRDS(results, file = 'Results.Rds')

#need to store cluster_num

#based on the output from the function
#write a second program - post processing program (the code needed to analyze the data from this function)
#handle each component of the results object to give you the effect among the different strata

str(Results)