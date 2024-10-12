This source code accompanies the report, detailing advice commissioned by the Australian Prime Minister and Cabinet, regarding the reopening of Australia in the late stages of the nation's COVID-19 vaccination campaign. 

This source code is a component of the model used to simulate infection importation events due to international travel. It simulates cohorts of arriving passengers moving through a specified quarantine system, and subsequently entering the community. This code does not simulate community transmission, it only accounts for transmission within the quarantine system. The output from this code (linelists of breach events), was used as input for subsequent modelling of community transmission scenarios. 

The code is implemented in Julia language. To run the code using the Juno IDE (Atom), 

set ...\families_version as the working directory, and enter the command: 

include("Quar_families_main_newVE_b015_2021_10_11.jl")

The above will simulate the scenarios defined in "scenario_table_VEi36_VEs67_2021_10_12a.csv"

to simulate the scenarios defined in "scenario_table_fam_followup_2021_10_20.csv"

enter the command: 

include("Quar_families_main_followup_b015_2021_10_20.jl")

to simulate the scenarios defined in "scenario_table_fam_followup_2021_10_21.csv"

enter the command: 

include("Quar_families_main_followup_b015_2021_10_21.jl").


This will create a directory for results, that will be populated with simulation outputs as the model iterates through the specified scenarios. 

If you receive the following error: 

"ERROR: LoadError: MethodError: no method matching initialise_agents! (" ... ") The applicable method may be too new"...

simply re-enter the include("...") command as listed above, the method should then initialise successfully. 

NOTE: the file Quar_ABM_utils.jl defines the variable tf (Float64), this is the number of days simulated, and determines the run time. By default tf = 2000000.0 days. This long time frame is set in order to simulate many breach events from which statistics can be drawn for initialisation of outbreak models. 

Further details of the model design, and explanation of each line-listed output value can be found here: 
https://arxiv.org/abs/2109.12799

NOTE: in the preprint article at https://arxiv.org/abs/2109.12799, the scenario specifics (i.e., testing frequency, duration of quarantine, and the quarantine pathway) are generic and do not apply to the modelling used in advice to the Australian government. The specifics of those scenarios are described in the released report, and are implemented in this source code as specified in the scenario files listed above.

NOTE: the scenario specified in the table: "scenario_table_fam_followup_2021_10_21.csv"
is for unvaccinated travellers. In the output of these simulations (defined in "Quar_families_main_followup_b015_2021_10_21.jl"), the line list contains an additional column "child", which takes values of 0 (adult) or 1 (child). This variable is not included in the previous versions because the identity of children was associated with their vaccination status (in the versions dated 2021_10_11 and 2021_10_08, all vaccinated individuals are adults, all unvaccinated individuals are children). 

