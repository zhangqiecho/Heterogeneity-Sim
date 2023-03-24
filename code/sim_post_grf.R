packages <- c("data.table","tidyverse","skimr","here")

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

readRDS(here("data","Results.Rds"))
#output from the first simulation run
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

saveRDS(results, file = here("data","Results.Rds"))

#need to store cluster_num

#based on the output from the function
#write a second program - post processing program (the code needed to analyze the data from this function)
#handle each component of the results object to give you the effect among the different strata

str(Results)