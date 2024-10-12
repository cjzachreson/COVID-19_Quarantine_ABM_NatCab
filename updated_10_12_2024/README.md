# COVID-19_Quarantine_ABM_NatCab
The directories contain four different versions of quarantine model simulation code. 
Each is derived from the model described in the arXiv preprint located at: 
https://arxiv.org/abs/2109.12799

The main differences are

base_version: same model, but implementing different testing schedules and isolating contacts of confirmed cases. 

families_version: same as base_version, but splitting families when a confirmed case is identified so that children are never isolated alone. 
In this version, children are not vaccinated 

no_intervention_version: same as base, but with minor modifications to make sure that all travellers enter the community immediately upon arrival 

RAT_version: reduced test sensitivity. 

More details for each version are contained in README files located in the separate directories. 
NOTE: the model is implemented in Julia language, and draws basic functionality from 
Agents.jl, a framework designed for agent-based modelling, see https://juliadynamics.github.io/Agents.jl/stable/

Update (October 12th, 2024): 
Fixed low-level error in within-host model to match calibration correctly. 
