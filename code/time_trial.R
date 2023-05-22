packages <- c("data.table","tidyverse","skimr","here")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}


a <- read_csv(here("data","run_time_1.csv"))
b <- read_csv(here("data","run_time_3.csv"))
c <- read_csv(here("data","run_time_9.csv"))

a
b
c

plot(rbind(a,b,c))

time_dat <- rbind(a,b,c)

# 1440 minutes per day
# predicted run time in days:
round((coef(lm(run_time ~ number_sims, data=time_dat))[2]*10000)/1440)

