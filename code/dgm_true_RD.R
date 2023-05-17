packages <- c("data.table","tidyverse","skimr","here","mvtnorm")

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

true_m0 <- true_m1 <- NULL
seq1 <- seq(-4,4,by=2)
seq2 <- seq(-6,6,by=2)
for(param1 in seq1){
  for(param2 in seq2){
    cluster_sim <- function(){
    
    # param1 = -4; param2 = -4
      
    param1 <- exp(param1)
    param2 <- exp(param2)
      
      n = 1e7
      p = 5
      set.seed(123)
      
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
      
      mu_y <- expit(-log((1/.2)-1) - log(param1)*pi_x - log(7.5)*pi_m + 
                      log(param2)*pi_x*pi_m + log(param1)*x+ log(7.5)*m - log(param2)*x*m)
      
      y <- rbinom(n, 1, mu_y)
      
      dd <- data.frame(y,m,c,x)
      names(dd) <- c("y","m",paste0("c",1:5),"x")
      
      return(dd)
    
    }
    
    sim_dat <- cluster_sim()
    
    mod_true <- lm(y ~ x + m + x*m, data = sim_dat)
    c_matrix1 <- c(0,1,0,0) # ATE among m = 0
    c_matrix2 <- c(0,1,0,1) # ATE among m = 1; NB: r puts interaction coefficient at end
    
    true_m0 <- rbind(true_m0,
                     data.frame(stratum = "m0",
                                estimate = coef(mod_true) %*% c_matrix1, 
                                param1 = param1, 
                                param2 = param2
                               )
                    )
    
    true_m1 <- rbind(true_m1,
                     data.frame(stratum = "m1",
                                estimate = coef(mod_true) %*% c_matrix2, 
                                param1 = param1, 
                                param2 = param2
                                )
                      )
    
  }
}

res <- rbind(true_m0,
             true_m1)

write_csv(res, file = here("data","true_RDs.csv"))