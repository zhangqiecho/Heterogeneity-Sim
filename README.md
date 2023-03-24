# Heterogeneity Sim

This repository contains code to conduct simulation studies on methods to estimate treatment heterogeneity. 

In this initial phase of the project, we are working on the `sim_run_grf.R` file. 

The purpose of this file in the initial phase is to:

- familarize ourselves with the generalized random forest function
- test the grf algorithm's ability to estimate CATEs in simple settings

The goal of the initial phase is the write programs that will enable us to deploy 
simulations on the (Rollins cluster)[https://ainaimi.github.io/resources/computing.html]

*Phase 1* will be completed when we have programs that can be deployed on the cluster
and that return figures and tables with results (bias, mse, CI coverage) from our 
preliminary (simple) dgms.

*Phase 2* will involve running the same programs with slightly more complex dgms and 
under different scenarios (e.g., sample size)

*Phase 3* will involve writing a separate set of programs that explore clustering 
algorithms for treatment heterogeneity (`sim_run_clust.R`)