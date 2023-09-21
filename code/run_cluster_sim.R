# run the function of cluster_sim using foreach with doParallel in HPC
start.time <- Sys.time()

packages <- c("data.table","tidyverse","SuperLearner","boot",
              "here","mvtnorm","foreach","doParallel","broom",
              "ggplot2","grf", "lmtest", "sandwich", "glmnet") # added a package

for (package in packages) {
  library(package, character.only=T)
}

func.packages <- c("data.table","tidyverse","SuperLearner","boot",
              "mvtnorm","broom", "grf", "lmtest", "sandwich")

phi_dat <- list(phi = seq(-2,2,by=2), 
                phi2 = seq(-2,2,by=2))

data_grid <- expand.grid(phi_dat, KEEP.OUT.ATTRS = FALSE)

cnum=5
number_sims=3
sample_size=1000
run_loop=1

source("cluster_sim_func.R")

#parallel::detectCores()

#n.cores <- parallel::detectCores() - 2 # the detected number of cores might not be active on HPC

n.cores = 20

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK" # qz: Fork doesn't work on windows
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()


par_res <- foreach(i = 1:nrow(data_grid), .packages = func.packages) %dopar% {
  
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


k <- 1:(nrow(data_grid)*number_sims)
blp_cf <- blp_drl <- sim_res <- NULL
for(i in 1:nrow(data_grid)){
  for(j in 1:number_sims){
    blp_cf[k]  <- list(par_res[[i]][[j]][[1]])
    blp_drl[k] <- list(blp_drl, par_res[[i]][[j]][[2]]) # qz: Q: why?
    sim_res[k] <- list(sim_res, par_res[[i]][[j]][[3]])
  }
}

blp_cf <- do.call(rbind,blp_cf)

blp_drl <- par_res[[1]][[1]][[2]]

sim_res <- par_res[[2]][[2]][[3]]

write_rds(par_res, 
          file = here("output",paste0("simulation_results_",number_sims,".rds"))
) # qz: changed the path

end.time <- Sys.time()

duration_time <- end.time - start.time

dur_dat <- data.frame(run_time = duration_time,
                      number_sims = number_sims)

write_csv(dur_dat, 
          file = here("output", 
                      paste0("run_time_",number_sims,run_loop,".csv"))) # qz: changed the path

