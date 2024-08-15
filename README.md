## GENERAL INFORMATION  
This README.txt file was updated on 2024-08-15.  

A. Manuscript associated with this archive:  
Mutz, J. & K.C. Abbott. 2024. Life history modulates the effects of inducible defenses on consumer-resource dynamics. American Naturalist.  
Abstract: Inducible defenses can affect the persistence, structure, and stability of consumer-resource systems. Theory shows that these effects depend on characteristics of the inducible defense, including timing, costs, efficacy, and sensitivity to consumer density. However, the expression and costs of inducible defenses often vary among life stages, which has not been captured in previous, unstructured models. To explore how inducible defenses expressed in stage-structured populations affect consumer-resource dynamics, we developed a model based on the biology of plant-herbivore interactions, with the plant (resource) population structured into juvenile and mature stages. We then investigated the joint effects of inducible defenses and resource life history (i.e., patterns of fecundity, maturation, and mortality) by simulating dynamics for plant populations occurring along a fast-slow pace-of-life continuum. In general, high inducible defense costs, or a slow pace-of-life coupled with high herbivore growth rates, promoted persistent cycles. However, these cycles fundamentally differed, with either the plant or the herbivore population peaking first. Additionally, plant population pace-of-life influenced the relative effects of stage-specific induction strength on equilibrium densities and the extent to which inducible defenses enabled persistence. Our work illustrates how life history modifies the population-level effects of trait-mediated interactions, with implications for conservation and pest management.  

B. Author details:  
Jessie Mutz, University of Tennessee, Knoxville (jmutz@utk.edu; corresponding author)  
Karen C. Abbott, Case Western Reserve University (kca27@case.edu)  

C. Funding sources:  
This work was funded by USDA-NIFA Award 2021-67034-35232.  


## ACCESS INFORMATION 
1. Licenses/restrictions placed on the code: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
2. The repository is available at: [http://doi.org/10.5281/zenodo.10501288] (http://doi.org/10.5281/zenodo.10501288).

## OVERVIEW OF SUPPLEMENTAL CODE AND SIMULATION OUTPUT  
This data repository consists of 3 code scripts (.R), 6 data files (.csv), and this README document.

### R scripts
1. **01_model-functions.R**  
The code in this script file: (i) defines the full ODE model as described in Table 2; (ii) defines function 'funcy' to simulate equilibrium dynamics across combinations of three parameters.

2. **02_simulations.R**  
The code in this script file: (i) samples random values of seedling maturation rate (g) and reproductive plant mortality rate (m) to use in simulations;  (ii) defines function 'sim_func' for running simulations across values of inducible defense responsiveness (tau), inducible defense cost (theta), and inducible defense effectiveness (epsilon); (iii) simulates the model to determine qualitative and quantitative behavior at equilibrium for various life history scenarios (NOT RUN -- output provided in R Workspace files 04-XX described below).

3. **03_figures.R**  
The code in this script file reproduces Figs 1-4 from the main text.

### CSV files
1. **04_tau-theta-ep27-qual.csv**
This file contains qualitative simulation output across values of inducible defense responsiveness (tau) and inducible defense costs (theta) for 500 plant life history scenarios (i.e., combination of seedling maturation rate [g] and reproductive plant mortality rate [m]).

2. **05_tau0.8-theta-ep-qual.csv**
This file contains qualitative simulation output across values of inducible defense costs (theta) and inducible defense efficacy (epsilon) for 500 plant life history scenarios (i.e., combination of seedling maturation rate [g] and reproductive plant mortality rate [m]).

3. **06_tau0.8-theta4-ep-g0.955-qual.csv**
This file contains qualitative simulation output across values of inducible defense efficacy (epsilon) at increased seedling maturation rate (g = 0.955) for 500 plant life history scenarios.

4. **07_tau-theta-ep27-r8-qual.csv**
This file contains qualitative simulation output across values of inducible defense responsiveness (tau) and inducible defense costs (theta) at high herbivore intrinsic growth rate (r = 8) for 500 plant life history scenarios (i.e., combination of seedling maturation rate [g] and reproductive plant mortality rate [m]).

The above data files (04-07) have the following structure:
(For detail on behavior diagnostics, e.g., 'pos', 'pt.eq', 'cyc1', see 'funcy' function in 01_model-functions.R)

| Column name | Description |   
| ----------- | ----------- | 
| rep | plant life history scenario (i.e., sampled combination of g and m values) |
| b, g, x, c, a, w, r, mu, m, theta, tau_s (tau_p), ep_s (ep_p) | parameter values associated with each simulation |
| S.pos, P.pos, I_s.pos, I_p.pos, H.pos | binary factor indicating whether the state variable has a positive value at equilibrium |
| S.neg, P.neg, I_s.neg, I_p.neg, H.neg | binary factor indicating whether negative values occur during simulation |
| S.pt.eq, P.pt.eq, I_s.pt.eq, I_p.pt.eq, H.pt.eq | binary factor indicating whether the state variable reaches a stable value at equilibrium |
| S.cyc1, P.cyc1, I_s.cyc1, I_p.cyc1, H.cyc1 | binary factor indicating cycling of the state variable at equilibrium |
| S.cyc2, P.cyc2, I_s.cyc2, I_p_cyc2, H.cyc2 | second indicator of cycling (coefficient of variation of densities through time) |
| S.NA, P.NA, I_s.NA, I_p.NA, I_s.NA, H.NA | binary factor indicating that NAs occur during simulation |
| qual | factor summarizing behavior of all state variables at equilibrium: 'extinct', 'pos.pt.eq' (positive point equilibrium), 'cycles', 'unk' (unknown -- output contains NAs) |

5. **08_tau-theta-ep-quan.csv**
This file contains quantitative simulation output across values of inducible defense responsiveness (tau), inducible defense costs (theta), and inducible defense efficacy (epsilon) for 4 representative plant life history scenarios (i.e., combination of seedling maturation rate [g] and reproductive plant mortality rate [m]).

6. **09_taus-taup-theta1.5-ep30-quan.csv**
This file contains quantitative simulation output across values of inducible defense responsiveness of seedlings (tau_s) and reproductive plants (tau_p) for 4 representative plant life history scenarios (i.e., combination of seedling maturation rate [g] and reproductive plant mortality rate [m]).

The above data files (08-09) have the following structure:

| Column name | Description |   
| ----------- | ----------- | 
| rep | plant life history scenario (i.e., sampled combination of g and m values) |
| b, g, x, c, a, w, r, mu, m, theta, tau_s (tau_p), ep_s (ep_p) | parameter values associated with each simulation |
| timestep | timestep of simulation |
| S, P, I_s, I_p, H | value of state variable at timestep t (note that I_s and I_p are totals, not means) |
| P.prop | proportion of plant biomass that is reproductive at timestep t  (i.e., P/(S + P)) |
| i_s_mean, i_p_mean | mean inducible defense level of seedling (i_s_mean) and reproductive plant (i_p_mean) biomass at timestep t (i.e., I_s/S or I_p/P) |
| H_SP | herbivore density (no. herbivores per plant biomass, i.e., H/(S + P)) at timestep t |
| seeds | no. seeds produced by reproductive plants at timestep t |

## SOFTWARE VERSIONS
Model simulations were run and figures produced in R v4.3.1 with the following package versions: deSolve v1.36, parallel v4.3.1, reshape2 v1.4.4, ggplot2 v3.4.2, gridExtra v2.3, ggtext v0.1.2, metR v0.14.0, doBy v4.6.17.
