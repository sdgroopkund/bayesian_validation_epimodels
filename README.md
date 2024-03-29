# Bayesian validation framework for dynamic epidemic models

This repository contains codes for replicating the analyses in the paper ‘Bayesian validation framework for dynamic epidemic models’.

These codes are broadly divided into two groups: (i) Codes for running the model; (ii) Codes for running the Bayesian validation framework.

# Codes for running the model. 

1. Model_n.R – Contains the function sim.mod.n that runs the HIV transmission model for a given set of parameters by solving the ODE equations using the C ODE solver rk4 

2. FitModel_n.R – Contains the function ‘fitmodel.n’ that fits the parameter 'tau' based on the community-specific and global parameters ('x' and 'theta'), fitted using least squares with R  optimization solver ‘optim’. See Section 3.1 in the draft for the definition of 'tau'. 

3. GetResults_simu.R – Contains the function ‘getresults.simu’ that runs the model on the fitted parameters and outputs HIV prevalence among other things for each community over time.

4. The folder named ‘C’ – Needed for running the model. Contain backend C codes to solve the ODEs. Should be in the same folder as the rest of these codes for the model to run.

# Codes for running the Bayesian validation framework 

1. Bayesian_validation_set1n3 – Codes to run Simulation Settings 1 (Null Model) and 3 (Faulty Model) in the paper

2. Bayesian_validation_set2 – Codes to run Simulation Settings 2 (Faulty Prior) in the paper

3. function_gen.R – Codes for different functions needed for Bayesian validation analysis

4. ggplot_plots.R – Codes to generate the plots given in the paper (discrepancy plots for Settings 1, 2 and 3 and the posterior distribution of 'theta 5' for Setting 2)

5. ptp_calc_mahalanobis.R – Codes to calculate the posterior tail probability using Mahalanobis distance for Settings 1, 2 and 3, as discussed in the paper.

The Codes for running the Bayesian validation framework as well as the R files to run the model (Model_n.R, FitModel_n.R, GetResults_simu.R) are hosted here. 
However, to run the simulation network one also needs to have codes in the 'C' folder (item 4 under Codes for running the model above), which contains codes to run the C ODE solver. To obtain these codes (to which we had no contribution to), one must reach out to the authors of Woods et al. (2018), the modeling team that built and calibrated this model.


# REFERENCES

1. Woods,  B.,  Rothery,  C.,  Anderson,  S.  J.,  Eaton,  J.  W.,  Revill,  P.,  Hallett,  T.  B.,  &  Claxton,  K.  (2018). Appraising the value of evidence generation activities: an HIV modelling study. BMJ  global  health, 3(6), e000488
