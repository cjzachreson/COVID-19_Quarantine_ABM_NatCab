# main quarantine ABM script

# turning detection off for benchmarking

#used seed = 10 on 2021 08 06

using DataFrames
using CSV
using Distributed


#global vaccinated_arrivals_index = [0, 1]
#global compliance_index = 1
#global efficacy_index = [1, 2]
#global R0_index = 0
#global pI_index = 0
#global R0s = 6.0#1.0:1.0:10.0
#global VEs = 0.0#0.0:0.1:0.9
#global pI_gs = [0.01]#:0.01:0.1
global R0 = 0.0
#global VE = 0.0
#global pI_g = 0.0
#global test_schedule_index = 1
#global compliance_rate_index = 1

#global scenario_filename = "scenario_table_C95_97_99_2021_09_23.csv"
global scenario_filename = "scenario_table_fam_followup_2021_10_20.csv"
scenario_list = DataFrame(CSV.File(scenario_filename, delim="\t"))
n_scenarios = size(scenario_list, 1)
global scenario_indices = 1:n_scenarios
global scenario_index = 0


function main()

    include("./Quar_ABM_utils.jl")
    include("./Agent_Quar.jl")
    include("./Quar_ABM_environment.jl")
    include("./Quar_ABM_disease_no_R0_b015.jl")
    include("./Quar_ABM_functions_2021_10_13.jl")

    date = "2021_10_20"

    output_dirname = pwd() * "/fam_split_followup_tf_2M/" * date *"/"
    if !isdir(output_dirname)
      #mkdir(output_dirname)
      mkpath(output_dirname)
    end

      model_main = AgentBasedModel(agent_type, Space)

      # initialise label for run:

      #tag = "R0_$R0"
      #tag = replace(tag, "." => "-")
      #label = label*tag*"_"

      #*******
      #declare default parameters
      #pVac_Travellers = 0.0
      #VEi = 0.0
      #VEs = 0.0
      #quar_duration = 0.0
      #test_schedule = [0.0]
      #compliance_rate = 0.0
      #hotel_quar_flag = false
      #home_quar_flag = false
      #sync_flag = false
      #*******

      tag,
      hotel_quar_flag,
      home_quar_flag,
      sync_flag,
      pVac_Travellers,
      VEi,
      VEs,
      quar_duration,
      test_schedule,
      compliance_rate =
      set_scenario_parameters!(scenario_index)

      #quar_duration = 14.1
      #test_schedule = [1.0, 5.0, 13.0]
      #quar_duration = 14.0

      if hotel_quar_flag
        quar_env = "Hotel"
      elseif home_quar_flag
        quar_env = "Home"
      else
        quar_env = "unspecified"
      end

      tag = "$tag" * "_" * "$quar_env"

      label = "Run_$scenario_index" * "_$tag" * "_"

      # set up agent properties:
      # test out Space by adding some agents:
      #set up some model parameters:
      n_workers = 20
      n_travellers = 100
      n_travellers_per_household = 4

      n_adults_per_household = 2
      n_children_per_household = 2 # must sum to n_travellers_per_household

      pVac_Workers = 1.0

      vaccinated_index_cases_flag = false

      # proportion of incomming travellers infected
      pI_0 = 0.01

      pI_Vac = pI_0 * (1 - VEs)

      FoI_fac_Hotel_TT_same = 1.0
      FoI_fac_Hotel_TT_dif = 0.01
      FoI_fac_Hotel_TW = 0.01
      FoI_fac_Hotel_WW = 0.1

      test_schedule_ext = test_schedule

      clinical_detection_flag = true

      quar_extension = 14.0
      iso_duration = 10.0
      test_schedule_ext_symptoms = [2.0, 8.0] # for contacts of those who express symptoms but cannot be isolated.
      iso_symptom_extension_flag = true

      # this makes sure groups moved
      home_ext_compliance_rate = 1.0

      #this means contacts of confirmed cases are taken to medi hotel for extension
      medi_hotel_flag = true # the rules on this will change in this version -
      # the medi hotel will apply to groups with confirmed cases that cannot be split

      if sync_flag
        t_burn = 0.0
      else
        t_burn = copy(t_burn_g)
      end
      # the params below are only used if sync_flag = true
      #***
      inter_cohort_interval = 7.0 # days between cohorts
      cohort_interval = quar_duration + inter_cohort_interval # may need to examine this choice...
      cohort_interval_clock = 0.0 # counter to separate cohorts
      #***

      test_report_delay = 1.0

      # select work roster:
      #roster = rosters[1] # everyone works 7 days per week
      #roster = rosters[2] # everyone works 5 days per week
      roster = rosters[3] # some people work 3 days, some work 5
      #roster = rosters[4] # some people work 1 day, some 3 days, and some work 5
      #TODO make sure correct number of workers take days off. not relevant to home quarantine.


      # generate config file for run:
      run_label = "$label" * date
      run_dirname = output_dirname * run_label
      println("recording output in file: $run_dirname")
      if !isdir(run_dirname)
        mkdir(run_dirname)
      end
      run_configfile_name = [run_dirname * "/$label" * "config.txt"]

      write_config_file(join(run_configfile_name),
          R0,
          tf,
          dt,
          t_burn,
          pI_0,
          pI_Vac,
          quar_duration,
          quar_extension,
          iso_duration,
          quar_env,
          sync_flag,
          inter_cohort_interval,
          cohort_interval,
          medi_hotel_flag,
          test_schedule,
          test_schedule_ext,
          test_report_delay,
          clinical_detection_flag,
          n_travellers_per_household,
          n_travellers,
          n_workers,
          roster,
          compliance_rate,
          pVac_Workers,
          pVac_Travellers,
          vaccinated_index_cases_flag,
          VEs,
          VEi,
          FoI_fac_Hotel_TT_same,
          FoI_fac_Hotel_TT_dif,
          FoI_fac_Hotel_TW,
          FoI_fac_Hotel_WW)


      #t_sens = []
      #i_check = []

      initialise_agents!(n_travellers, n_workers, model_main, environments)

      # check work schedule
      n_min = 5 # need a schedule with at least 5 workers present each day
      global work_schedule_OK = false;
      global trys = 0
      while !work_schedule_OK
        initialise_work_schedule!(model_main, environments, roster)
        global trys += 1 # julia while loop weirdness....
        global work_schedule_OK = check_work_schedule(model_main, n_min)
        if trys > 1000
          println("having trouble assigning minimum number of workers each day - adjust roster or number of workers")
          return
        end
      end

      # assign household IDs to travellers

      assign_households!(model_main,
                         n_travellers_per_household,
                         households,
                         n_children_per_household,
                         n_adults_per_household )
      # removal from memory occurs on a household basis. household groups are removed and replaced together.
      # note: removal is distinct from discharge - if members of a household are in iso,
      # the rest of the family can enter the community while waiting for them.

      # vaccinate agents - can swap this around with infection algorithm if desired.
      # this is now modified - children are not vaccinated
      vaccinate_all!(model_main, pVac_Workers, pVac_Travellers)


      # infect arrivals -
      # anyone in the arrivals node can become infected - remember to move them out
      # of arrivals node after evaluating infection status...
      infect_arrivals!(model_main, pI_0, pI_Vac)

      if vaccinated_index_cases_flag
        vaccinate_infected_arrivals!(model_main)
      end


      assign_compliance_rate_to_arrivals!(model_main, compliance_rate)

      # move arrivals into quarantine
      if hotel_quar_flag
        move_arrivals_to_hotel_quarantine!(model_main, test_schedule)
      elseif home_quar_flag
        move_arrivals_to_home_quarantine!(model_main, test_schedule)
      end
      #move_arrivals_to_hotel_quarantine!(model, test_schedule)


      # begin iterating through time.
      # infect a worker and see if testing works:
      #infect_agent!(model[1])

      global n_travellers_discharged = 0.0
      global n_infected_arrivals = 0.0


      #build weekly data vectors
      travellers_discharged_weekly = []
      infectious_travellers_dischaged_weekly = []
      secondary_cases_weekly = []
      exposure_days_T_weekly = []
      exposure_days_T_uv_weekly = []
      exposure_days_T_vac_weekly = []
      exposure_days_W_weekly = []
      exposure_days_W_uv_weekly = []
      exposure_days_W_vac_weekly = []
      infected_arrivals_uv_weekly = []
      infected_arrivals_vac_weekly = []
      PCR_detections_weekly = []
      clinical_detections_weekly = []

      travellers_discharged_this_week = 0.0
      infectious_travellers_dischaged_this_week = 0.0
      secondary_cases_this_week = 0.0
      exposure_days_T_this_week = 0.0
      exposure_days_T_uv_this_week = 0.0
      exposure_days_T_vac_this_week = 0.0
      exposure_days_W_this_week = 0.0
      exposure_days_W_uv_this_week = 0.0
      exposure_days_W_vac_this_week = 0.0
      infected_arrivals_uv_this_week = 0.0
      infected_arrivals_vac_this_week = 0.0
      PCR_detections_this_week = 0.0
      clinical_detections_this_week = 0.0


      time_increment = copy(dt)

      #initialise list of new households (necessary if sync_flag is on)
      n_households_to_add = [0]
      n_agents_to_add = [0]

      for t in 1.0:time_increment:tf

          if mod(t, 1000) == 0
            println(t)
          end

          if mod(t, 1) == 0
            #initialise output summary tallies for the first week:

            # record weekly summary data
            if mod(t, 7) == 0
              if t > t_burn
                #println("recording at time $t, $(mod(t, 7))")
                push!(travellers_discharged_weekly, travellers_discharged_this_week)
                push!(infectious_travellers_dischaged_weekly, infectious_travellers_dischaged_this_week)
                push!(secondary_cases_weekly, secondary_cases_this_week)
                push!(exposure_days_T_weekly, exposure_days_T_this_week)
                push!(exposure_days_T_uv_weekly, exposure_days_T_uv_this_week)
                push!(exposure_days_T_vac_weekly, exposure_days_T_vac_this_week)
                push!(exposure_days_W_weekly, exposure_days_W_this_week)
                push!(exposure_days_W_uv_weekly, exposure_days_W_uv_this_week)
                push!(exposure_days_W_vac_weekly, exposure_days_W_vac_this_week)
                push!(infected_arrivals_uv_weekly, infected_arrivals_uv_this_week)
                push!(infected_arrivals_vac_weekly, infected_arrivals_vac_this_week)
                push!(PCR_detections_weekly, PCR_detections_this_week)
                push!(clinical_detections_weekly, clinical_detections_this_week)

                travellers_discharged_this_week = 0.0
                infectious_travellers_dischaged_this_week = 0.0
                secondary_cases_this_week = 0.0
                exposure_days_T_this_week = 0.0
                exposure_days_T_uv_this_week = 0.0
                exposure_days_T_vac_this_week = 0.0
                exposure_days_W_this_week = 0.0
                exposure_days_W_uv_this_week = 0.0
                exposure_days_W_vac_this_week = 0.0
                infected_arrivals_uv_this_week = 0.0
                infected_arrivals_vac_this_week = 0.0
                PCR_detections_this_week = 0.0
                clinical_detections_this_week = 0.0
              end
            end
          end


          n_workers_present = convert(Float64, n_workers)


        #  evaluate_traveller_tests!(model_main)
          # check to see if it's a new day:
          if mod(t, 1) == 0
            #println(length(households))
            # 1) the workforce for the day is determined,
            local day_of_week = mod(t, 7)
            n_workers_present = assign_workforce!(model_main, day_of_week)
            #2) the workforce present is tested
            test_workforce!(model_main, test_report_delay)
            # 3) test results for workers previously tested are returned and evaluated
            evaluate_workforce_tests!(model_main)
            # 4) positive tests or onset of symptoms in workers cause them to be 'discharged'
            # 5) a discharged worker is replaced with a new, susceptible worker with the same schedule
            discharge_workers!(model_main, t, pVac_Workers, p_asymp)

            #6) test results are returned for travellers previously tested
            evaluate_traveller_tests!(model_main)
            # 7) travellers who test positive or recently became symptomatic are put into isolation
            # 8) contacts of travellers who test positive are put into extension (if they're not already there)

            PCR_detections_t,
            clinical_detections_t =
            isolate_travellers_family!(model_main,
                                       iso_duration, quar_extension,
                                       test_schedule_ext,
                                       medi_hotel_flag,
                                       home_ext_compliance_rate,
                                       test_schedule_ext_symptoms)



            PCR_detections_this_week += PCR_detections_t
            clinical_detections_this_week += clinical_detections_t
            # 9) all travellers who are finished with quarantine are discharged into the community

            # synchronise re-initialisation of households for cohorting.
            if sync_flag
                cohort_interval_clock += 1.0
            else
                n_households_to_add = [0]
                n_agents_to_add = [0]
            end

            n_discharged, n_discharged_infected = discharge_travellers!(model_main, t, n_households_to_add, n_agents_to_add)
            # 10) once all members of a household are discharged, a new household is added to the arrivals node
            travellers_discharged_this_week += n_discharged
            infectious_travellers_dischaged_this_week += n_discharged_infected

            if sync_flag
                if cohort_interval_clock >= cohort_interval
                    # new arrivals are initialised (vaccine and infection status is determined)
                    #println("adding new cohort of size $n_agents_to_add at time $t")
                    #println("interval clock: $cohort_interval_clock")
                    n_infected_arrivals_t,
                    n_infected_arrivals_uv_t,
                    n_infected_arrivals_vac_t =
                    initialise_arrivals!(model_main,
                                        pVac_Travellers, pI_Vac,
                                        pI_0, n_households_to_add[1],
                                        n_agents_to_add[1], test_schedule,
                                        n_children_per_household,
                                        n_adults_per_household)

                    infected_arrivals_uv_this_week += n_infected_arrivals_uv_t
                    infected_arrivals_vac_this_week += n_infected_arrivals_vac_t

                    # zero the sync lists
                    cohort_interval_clock = 0.0
                    n_households_to_add = [0]
                    n_agents_to_add = [0]

                    #println("zeroing interval clock at time $t")
                    #println("there are: $(nagents(model_main)) agents" )

                end
            else
                n_infected_arrivals_t,
                n_infected_arrivals_uv_t,
                n_infected_arrivals_vac_t =
                initialise_arrivals!(model_main,
                                     pVac_Travellers, pI_Vac,
                                     pI_0, n_households_to_add[1],
                                     n_agents_to_add[1], test_schedule,
                                     n_children_per_household,
                                     n_adults_per_household)

                infected_arrivals_uv_this_week += n_infected_arrivals_uv_t
                infected_arrivals_vac_this_week += n_infected_arrivals_vac_t
            end

            # makes sure all index cases are vaccinated.
            if vaccinated_index_cases_flag
              vaccinate_infected_arrivals!(model_main)
            end
            assign_compliance_rate_to_arrivals!(model_main, compliance_rate)

            # the below functions can be modified to move conditional on agent properties
            if hotel_quar_flag
              move_arrivals_to_hotel_quarantine!(model_main, test_schedule)
            elseif home_quar_flag
              move_arrivals_to_home_quarantine!(model_main, test_schedule)
            end

            #println(nagents(model_main))

            # once all individuals in a household are in the community and recovered,
            # they are removed from the model. NOTE: the delete! function can be used delete
            # the specified households from the households dict. and the kill_agent! function
            # can be used to remove the agents from the model.

            record_and_delete!(model_main)

            test_travellers!(model_main, test_report_delay)

            #println(t)
            #check_test_sensitivity!(model, t_sens, i_check)

            evaluate_compliance!(model_main)

          end

          # TODO: compute transmission for different contact management policies

          secondary_cases_t = compute_transmission!(model_main,
                                                    VEi,
                                                    VEs,
                                                    FoI_fac_Hotel_TT_same,
                                                    FoI_fac_Hotel_TT_dif,
                                                    FoI_fac_Hotel_TW,
                                                    FoI_fac_Hotel_WW,
                                                    n_workers_present,
                                                    time_increment)

          secondary_cases_this_week += secondary_cases_t

          # step clocks and update states:
          edays_W_t,
          edays_uv_W_t,
          edays_vac_W_t = step_workers!(model_main, time_increment)

          exposure_days_W_this_week += edays_W_t
          exposure_days_W_uv_this_week += edays_uv_W_t
          exposure_days_W_vac_this_week += edays_vac_W_t



          edays_T_t,
          edays_uv_T_t,
          edays_vac_T_t = step_travellers!(model_main,
                                           quar_duration,
                                           iso_duration,
                                           time_increment)

         exposure_days_T_this_week += edays_T_t
         exposure_days_T_uv_this_week += edays_uv_T_t
         exposure_days_T_vac_this_week += edays_vac_T_t


        if clinical_detection_flag
          assess_symptoms!(model_main, iso_duration, iso_symptom_extension_flag)
        end

        #check_agents!(model_main, t)


      end

      output_filename_T = run_dirname * "/Traveller_breach.csv"
      output_filename_W = run_dirname * "/Worker_breach.csv"

      outname_full_T = run_dirname * "/Travellers_full.csv"

      breach_T = linelist_T[linelist_T.exposure_days .!= 0, :]
      breach_W = linelist_W

      CSV.write(output_filename_T, breach_T)
      CSV.write(output_filename_W, breach_W)

      CSV.write(outname_full_T, linelist_T)

      #println("discharged $n_travellers_discharged travellers" )
      #println("infected arrivals: $n_infected_arrivals" )

      #output_filename_tsens = output_dirname * "/tsens.csv"
      #CSV.write(output_filename_tsens, DataFrame(temp_ = t_sens))

      # weekly averages:
      if sync_flag
          w1_stst = 1
      else
        w1_stst = convert(Int64, ceil(t_burn / 7))
      end

      travellers_discharged_pw = mean(travellers_discharged_weekly[w1_stst:end])
      I_travellers_discharged_pw = mean(infectious_travellers_dischaged_weekly[w1_stst:end])
      secondary_cases_pw = mean(secondary_cases_weekly[w1_stst:end])
      edays_pw_T = mean(exposure_days_T_weekly[w1_stst:end])
      edays_pw_T_uv = mean(exposure_days_T_uv_weekly[w1_stst:end])
      edays_pw_T_vac = mean(exposure_days_T_vac_weekly[w1_stst:end])
      edays_pw_W = mean(exposure_days_W_weekly[w1_stst:end])
      edays_pw_W_uv = mean(exposure_days_W_uv_weekly[w1_stst:end])
      edays_pw_W_vac = mean(exposure_days_W_vac_weekly[w1_stst:end])
      infected_arrivals_pw_uvac = mean(infected_arrivals_uv_weekly[w1_stst:end])
      infected_arrivals_pw_vac = mean(infected_arrivals_vac_weekly[w1_stst:end])
      detections_pw_tests = mean(PCR_detections_weekly[w1_stst:end])
      detections_pw_clin = mean(clinical_detections_weekly[w1_stst:end])

      edays_per_inf_arrival_T = edays_pw_T / (infected_arrivals_pw_uvac + infected_arrivals_pw_vac)
      edays_per_inf_arrival_W = edays_pw_W / (infected_arrivals_pw_uvac + infected_arrivals_pw_vac)

      # generate results summary file for run:

      results_summary_filename = join([run_dirname * "/$label" * "results_summary.txt"])

      write_output_summary(
        results_summary_filename,
        edays_pw_T,
        edays_pw_T_uv,
        edays_pw_T_vac,
        edays_pw_W,
        edays_pw_W_uv,
        edays_pw_W_vac,
        edays_per_inf_arrival_T,
        edays_per_inf_arrival_W,
        infected_arrivals_pw_uvac,
        infected_arrivals_pw_vac,
        detections_pw_tests,
        detections_pw_clin,
        travellers_discharged_pw,
        I_travellers_discharged_pw,
        secondary_cases_pw)

end

for s_i in 1:n_scenarios
      global R0 = 6.0
      global scenario_index = scenario_indices[s_i]

      main()
end





# every day, the following events happen:
# 1) the workforce for the day is determined,
# 2) the workforce present is tested
# 3) test results for workers previously tested are returned and evaluated
# 4) positive tests or onset of symptoms in workers cause them to be 'discharged'
# 5) a discharged worker is replaced with a new, susceptible worker with the same schedule

# 6) test results are returned for travellers previously tested
# 7) travellers who test positive or recently became symptomatic are put into isolation
# 8) contacts of travellers who test positive are put into extension (if they're not already there)
# 9) all travellers who are finished with quarantine are discharged into the community
# 10) once all members of a household are discharged, a new household is added in the 'arrivals' node and initialised
    # once all members of a discharged household are recovered, they are tabulated into the line list and removed from memory
# 12) infection and vaccination status of new arrivals is evaluated and they are moved into quarantine
      # optionally, this can be handled by the 'step travellers' function
# 13) all travellers scheduled for testing are tested (results returned after specified delay)
      # note - test sensitivity is evaluated based on infection clock and individual-level parameters

# every timestep, the following events happen:
# 12) iterate through the occupied nodes of the graph and compute pair-wise transmission between contacts
# 13) step all clocks
# TODO decide whether symptomatic agents should be isolated immediately, or on the following day
