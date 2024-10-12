# utilities

using Random
using DataFrames
using CSV

# random number generator for system:
#seed = 10 # used on and before 2021 09 21
seed = 20 # altered for overnight run 2021 09 22
rng = MersenneTwister(seed);

global dt = 0.1
global tf = 100000.0
#global tf = 2000.0
#global tf = 1000000.0
#global tf = 2000000.0

global t_burn_g = 1000.0

global linelist_W = DataFrame(exposure_days = [0.0],
                       t_latent = [0.0],
                       t_incubation = [0.0],
                       t_post_incubation = [0.0],
                       tested_positive = [0],
                       expressed_symptoms = [0],
                       time_discharged = [0.0],
                       symptomatic = [0],
                       vaccinated = [0],
                       FoI_max = [0.0],
                       FoI_community = [0.0]
                      )

global linelist_W_entry = Dict(:exposure_days => 0.0,
                       :t_latent => 0.0,
                       :t_incubation => 0.0,
                       :t_post_incubation => 0.0,
                       :tested_positive => 0,
                       :expressed_symptoms => 0,
                       :time_discharged => 0.0,
                       :symptomatic => 0,
                       :vaccinated => 0,
                       :FoI_max => 0.0,
                       :FoI_community => 0.0
                      )

global linelist_T = DataFrame(exposure_days = [0.0],
                       days_in_quar = [0.0],
                       days_in_extended_quar = [0.0],
                       days_in_iso = [0.0],
                       t_latent = [0.0],
                       t_incubation = [0.0],
                       t_post_incubation = [0.0],
                       time_discharged = [0.0],
                       index_case = [0],
                       symptomatic = [0],
                       vaccinated = [0],
                       compliance = [0.0],
                       FoI_max = [0.0],
                       FoI_community = [0.0]
                      )

global linelist_T_entry = Dict(:exposure_days => 0.0,
                        :days_in_quar => 0.0,
                        :days_in_extended_quar => 0.0,
                        :days_in_iso => 0.0,
                        :t_latent => 0.0,
                        :t_incubation => 0.0,
                        :t_post_incubation => 0.0,
                        :time_discharged => 0.0,
                        :index_case => 0,
                        :symptomatic => 0,
                        :vaccinated => 0,
                        :compliance => 0.0,
                        :FoI_max => 0.0,
                        :FoI_community => 0.0
                      )

global scenarios = DataFrame(CSV.File(scenario_filename, delim="\t"))

function write_config_file(
  output_filename::String,
  R0::Float64,
  t_final::Float64,
  delta_t::Float64,
  t_relax::Float64,
  imporation_rate::Float64,
  vac_importation_rate::Float64,
  quar_duration::Float64,
  ext_duration::Float64,
  iso_duration::Float64,
  quar_environment::String,
  sync_flag::Bool,
  inter_cohort_interval::Float64,
  cohort_interval::Float64,
  medi_hotel_flag::Bool,
  test_days::Vector{Float64},
  ext_test_days::Vector{Float64},
  test_report_delay::Float64,
  clinical_detection::Bool,
  groups_size::Int64,
  n_travellers::Int64,
  n_workers::Int64,
  work_roster::AbstractVector,
  compliance_prob::Float64,
  pVac_Workers::Float64,
  pVac_Travellers::Float64,
  Vaccinated_index_cases::Bool,
  Vac_efficacy_infection::Float64,
  Vac_efficacy_transmission::Float64,
  transmission_Fac_TTsame::Float64,
  transmissionFac_TTdiff::Float64,
  transmission_Fac_TW::Float64,
  transmission_Fac_WW::Float64)

  open(output_filename, "w") do io
    write(io, join(["configuration input parameters: \n
             ******* \n
             Reproductive Ratio (homogeneous calib.): $R0  \n
             simulation time: $t_final  \n
             time step: $delta_t  \n
             burn-in time (arb): $t_relax \n
             baseline importation rate p(I | unvac., arrival): $imporation_rate \n
             vaccinated importation rate p(I | vac., arrival): $vac_importation_rate \n
             quarantine duration: $quar_duration \n
             contact extension duration: $ext_duration \n
             case isolation duration (extended from symptom onset): $iso_duration \n
             quarantine environment: $quar_environment \n
             synchronised cohorts: $sync_flag \n
             clearning period between cohorts (if sync): $inter_cohort_interval \n
             interval between cohorts (if sync): $cohort_interval \n
             contacts transferred to MediHotel: $medi_hotel_flag \n
             test days: $test_days \n
             extension test days: $ext_test_days \n
             test report delay: $test_report_delay \n
             clinical detection (via symptoms): $clinical_detection \n
             group size (close contacts): $groups_size \n
             n travellers: $n_travellers \n
             n workers: $n_workers \n
             work roster: $work_roster \n
             compliance probability: $compliance_prob \n
             proportion of workers vaccinated: $pVac_Workers \n
             proportion of travellers vaccinated: $pVac_Travellers \n
             all index cases vaccinated: $Vaccinated_index_cases \n
             vaccine efficacy against infection: $Vac_efficacy_infection \n
             vaccine efficacy against transmission: $Vac_efficacy_transmission \n
             Hotel transmission factor (Traveller, Traveller, same group): $transmission_Fac_TTsame \n
             Hotel transmission factor (Traveller, Traveller, different group): $transmissionFac_TTdiff \n
             Hotel transmission factor (Traveller, Worker): $transmission_Fac_TW \n
             Hotel transmission factor (Worker, Worker): $transmission_Fac_WW \n
             *******"]))
           end
end

function write_output_summary(
    output_filename::String,
    edays_pw_T::Float64,
    edays_pw_T_uv::Float64,
    edays_pw_T_vac::Float64,
    edays_pw_W::Float64,
    edays_pw_W_uv::Float64,
    edays_pw_W_vac::Float64,
    edays_per_inf_arrival_T::Float64,
    edays_per_inf_arrival_W::Float64,
    infected_arrivals_pw_uvac::Float64,
    infected_arrivals_pw_vac::Float64,
    detections_pw_tests::Float64,
    detections_pw_clin::Float64,
    travellers_discharged_pw::Float64,
    I_travellers_discharged_pw::Float64,
    secondary_cases_pw::Float64)

  open(output_filename, "w") do io
    write(io, join(["simulation output summary: \n
             ******* \n
             exposure days per week (Travellers): $edays_pw_T
             exposure days per week (Travellers, unvaccinated): $edays_pw_T_uv
             exposure days per week (Travellers, vaccinated): $edays_pw_T_vac \n
             exposure days per week (Workers): $edays_pw_W \n
             exposure days per week (Workers, unvaccinated): $edays_pw_W_uv \n
             exposure days per week (Workers, vaccinated): $edays_pw_W_vac \n
             exposure days per infected arrival (Travellers): $edays_per_inf_arrival_T \n
             exposure days per infected arrival (Workers): $edays_per_inf_arrival_W \n
             infected arrivals per week (unvaccinated): $infected_arrivals_pw_uvac \n
             infected arrivals per week (vaccinated): $infected_arrivals_pw_vac \n
             detected cases per week (PCR test): $detections_pw_tests \n
             detected cases per week (symptom onset): $detections_pw_clin \n
             travellers discharged per week: $travellers_discharged_pw \n
             infected travellers discharged per week: $I_travellers_discharged_pw \n
             secondary cases per week: $secondary_cases_pw \n
             *******"]))
           end
end

function set_scenario_parameters!(i)::Tuple{String, Bool, Bool, Bool, Float64, Float64, Float64, Float64, Vector{Float64}, Float64}

    quar_duration = convert(Float64, scenarios.quar_duration[i])
    tag = "D" * "$quar_duration"

    test_schedule = parse.(Float64, split(scenarios.test_schedule[i], ", "))
    tag = "$tag" * "_T" * "$(size(test_schedule, 1))"

    home_quar_flag = scenarios.home_quar_flag[i]
    hotel_quar_flag = scenarios.hotel_quar_flag[i]

    sync_flag = scenarios.sync_flag[i]
    if sync_flag
      tag = "$tag" * "_S"
    end

    pVac_Travellers = convert(Float64, scenarios.pVac_Travellers[i])
    tag = "$tag" * "_pV" * "$pVac_Travellers"

    if pVac_Travellers != 0.0
      VEi = scenarios.VEi[i]
      tag = "$tag" * "_Vi" * "$VEi"

      VEs = scenarios.VEs[i]
      tag = "$tag" * "_Vs" * "$VEs"
    else
      VEi = 0.0
      VEs = 0.0
    end

    compliance_rate = convert(Float64, scenarios.compliance_rate[i])
    tag = "$tag" * "_C" * "$compliance_rate"

    tag = replace("$tag", "." => "")
    tag = replace("$tag", "0" => "")

    println(tag)

    return tag, hotel_quar_flag, home_quar_flag, sync_flag, pVac_Travellers, VEi, VEs, quar_duration, test_schedule, compliance_rate

end
