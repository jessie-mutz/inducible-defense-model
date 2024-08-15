## GENERAL INFORMATION  
This README.txt file was updated on 2024-01-13.  
A. Manuscript associated with this archive: Life history modulates the effects of inducible defenses on consumer-resource dynamics.  
All other general information omitted for double-blind review.

## ACCESS INFORMATION 
1. Licenses/restrictions placed on the code: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
2. Recommended citation for this code archive: 

## OVERVIEW OF SUPPLEMENTAL CODE AND SIMULATION OUTPUT  
This data repository consists of 6 code scripts, 4 R workspaces, and this README document.

### R scripts
1. **01_model-functions.R**  
The code in this script file: (i) defines the full ODE model as described in Table S1; (ii) defines functions 'run_qual' and 'qual_out' to simulate qualitative behavior across combinations of two parameters; (iii) defines function 'funcy' to simulate equilibrium dynamics across combinations of three parameters.

2. **02_qualitative-simulations.R**  
The code in this script file: (i) defines function 'qual_set' for running qualitative behavior simulations for tau x epsilon, tau x b, tau x g, and tau x m; (ii) simulates the model for four plant life history scenarios (output provided in R Workspace files 07-10 described below).

3. **03_qualitative-figures.R**  
The code in this script file: (i) uses qualitative behavior simulation output to define regions of parameter space (tau x epsilon, tau x theta, tau x r, or tau x m) where the system goes extinct, goes to a stable point equilibrium, or has persistent cycles; (ii) reproduces Figs 1 & S1.

4. **04_cycling-simulations-and-figures.R**  
The code in this script file: (i) simulates equilibrium dynamics reported in Figs 2-3 & S2-5; (ii) reproduces Figs 2-3 & S2-5.

5. **05_quantitative-simulations-and-figures.R**  
The code in this script file: (i) simulates equilibrium dynamics for four plant life history scenarios, across values of tau, theta, and epsilon; (ii) reproduces Figs 4 & S6-7.

6. **06_uneven-tau-simulations-and-figures.R**  
The code in this script: (i) simulates equilibrium dynamics for four plant life history scenarios, across values of tau_S and tau_P; (ii) reproduces Figs 5 & S8-9.


### R Workspaces
Each of the following workspaces contains 4 lists: 'tau_ep_out', 'tau_b_out', 'tau_g_out', and 'tau_m_out'. 

| Object | Description |  
| ------ | ----------- |  
| tau_ep_out | qualitative behavior at equilibrium across values of tau (induction responsiveness) and epsilon (induction effectiveness) |
| tau_b_out | qualitative behavior at equilibrium across values of tau (induction responsiveness) and b (plant fecundity) |
| tau_g_out | qualitative behavior at equilibrium across values of tau (induction responsiveness) and g (seedling maturation rate) |
| tau_m_out | qualitative behavior at equilibrium across values of tau (induction responsiveness) and m (mature plant mortality rate) |

Each list has named elements: 'baseline', 'theta = 4', 'c = 5', and 'r = 8', with the structure:  
(For detail on behavior diagnostics, e.g., 'pos', 'pt.eq', 'cyc1', see 'run_qual' function in 01_model-functions.R)

| Column name | Description |   
| ----------- | ----------- | 
| b, theta, g, x, tau_s, tau_p, c, a, ep_s, ep_p, w, r, mu, m | parameter values associated with each simulation |
| S.pos, P.pos, I_s.pos, I_p.pos, H.pos | binary factor indicating whether the state variable has a positive value at equilibrium |
| S.pt.eq, P.pt.eq, I_s.pt.eq, I_p.pt.eq, H.pt.eq | binary factor indicating whether the state variable reaches a stable value at equilibrium |
| S.cyc1, P.cyc1, I_s.cyc1, I_p.cyc1, H.cyc1 | binary factor indicating cycling of the state variable at equilibrium |
| S.cyc2, P.cyc2, I_s.cyc2, I_p_cyc2, H.cyc2 | second indicator of cycling (coefficient of variation of densities through time) |
| S.neg, P.neg, I_s.neg, I_p.neg, H.neg | binary factor indicating that negative values occur during simulation |
| S.NA, P.NA, I_s.NA, I_p.NA, I_s.NA, H.NA | binary factor indicating that NAs occur during simulation |
| sum.pos | sum of 'pos' indicators across the five state variables |
| sum.pt.eq | sum of 'pt.eq' indicators across the five state variables |
| sum.cyc1 | sum of 'cyc1' indicators across the five state variables |
| sum.NA | sum of 'NA' indicators across the five state variables |
| sum.pos.pt.eq | sum of 'sum.pos' and 'sum.pt.eq', used to determine whether all five state variables reach a positive point equilibrium |
| qual | factor summarizing behavior of all state variables at equilibrium: 'extinct', 'pos.pt.eq' (positive point equilibrium), 'cycles', 'unk' (unknown -- output contains NAs) |

1. **07_qualitative-output-lowm-lowg.RData**  
Simulation output for first plant life history scenario: b = 6, g = 0.15, m = 0.05 (see Methods in main text & Table S2).

3. **08_qualitative-output-highm-lowg.RData**  
Simulation output for second plant life history scenario: b = 11, g = 0.15, m = 0.8 (see Methods in main text & Table S2).

4. **09_qualitative-output-lowm-modg.RData**  
Simulation output for third plant life history scenario: b = 2.43, g = 0.5, m = 0.05 (see Methods in main text & Table S2). 

5. **10_qualitative-output-highm-modg.RData**  
Simulation output for second plant life history scenario: b = 4.445, g = 0.5, m = 0.8 (see Methods in main text & Table S2). 


## SOFTWARE VERSIONS
Model simulations were run and figures produced in R v4.3.1 with the following package versions: deSolve v1.36, parallel v4.3.1, reshape2 v1.4.4, ggplot2 v3.4.2, gridExtra v2.3, ggtext v0.1.2, metR v0.14.0, doBy v4.6.17.
