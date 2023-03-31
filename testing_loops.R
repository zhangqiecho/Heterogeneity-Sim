packages <- c("data.table","tidyverse","skimr","here","rbenchmark")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

benchmark("for_list" = {
          my_list <- list()
          my_list
          for(i in 1:2000) {
            output <- rep(i, 5)
            my_list[[i]] <- output
          }
        },
        "rbind_list" = {
          my_list <- NULL
          for(i in 1:2000) {
            output <- rep(i, 5)
            my_list <- rbind(my_list,output)
          }
        },
        replications = 5000
  )
