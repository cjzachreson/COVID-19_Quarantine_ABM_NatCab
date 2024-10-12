# definition of Agent structure


using Agents


# define household
mutable struct household_type
  id::Int
  members::Array{Int}
  is_split::Bool
end


# define agents:
mutable struct agent_type <: AbstractAgent
  id::Int
  pos::Int
  # define extra attributes
  # category of agent:
  worker::Bool
  traveller::Bool
  # traveller states
  is_child::Bool
  is_adult::Bool

  #worker-specific properties
  workdays::Array{Int64}
  present_in_workplace::Bool
  # traveller-specific properties:
  hh_id::Int64
  p_compliance::Float64 # probability of being compliant on a given day
  compliant::Bool
  quarantine_extended::Bool
  isolation_extended::Bool
  # some clocks
  time_in_quar::Float64
  time_in_extended_quar::Float64
  time_remaining_in_ext_quar::Float64 # new in this version - because ext. quar can be extended further.
  extension_group::Int64 # a number that gives a group id for split families.
  time_in_isolation::Float64
  time_remaining_in_isolation::Float64
  #general properties
  vaccinated::Bool
  infected::Bool
  latent::Bool
  incubating::Bool
  recovering::Bool
  recovered::Bool
  asymptomatic::Bool
  expressing_symptoms::Bool
  newly_symptomatic::Bool
  # some clocks
  time_infected::Float64
  time_incubating::Float64
  time_recovering::Float64
  time_infectious_in_community::Float64
  test_sensitivity::Float64
  b1::Float64
  b2::Float64
  b3::Float64
  breakpoint::Float64
  beta_t::Float64
  beta_max::Float64
  beta_min::Float64
  k_inc::Float64
  k_r::Float64

  incubation_period::Float64
  recovery_period::Float64

  test_positive::Bool
  awaiting_test_result::Bool
  time_to_test_result::Float64

  discharged::Bool
  detected::Bool
  detected_clinical::Bool
  detected_test::Bool
  already_detected::Bool
  recorded::Bool

  FoI_in_community::Float64

  time_discharged::Float64

  time_until_test::Array{Float64}

  index_case::Bool

  agent_type(id, pos) = new(id, pos,
                            false, #worker::Bool
                            false, #traveller::Bool
                            false, #is_child::Bool
                            false, #is_adult::Bool
                            [0], #work days
                            false, # present_in_workplace::Bool
                            0, #hh_id::Int64
                            1.0,  #p_compliance::Float64
                            true, #compliant::Bool
                            false, #quarantine_extended::Bool
                            false, #isolation_extended::Bool
                            0.0, #time_in_quar::Float64
                            0.0, #time_in_extended_quar::Float64
                            0.0, #time_remaining_in_ext_quar::Float64
                            0,   #extension_group::Int64
                            0.0, #time_in_isolation::Float64
                            0.0,  #time_remaining_in_isolation::Float64
                            false, #vaccinated::Bool
                            false, #infected::Bool
                            false, #latent::Bool
                            false, #incubating::Bool
                            false, #recovering::Bool
                            false, #recovered::Bool
                            false, #asymptomatic::Bool
                            false, #expressing_symptoms::Bool
                            false, #newly_symptomatic::Bool
                            0.0, #time_infected::Float64
                            0.0, #time_incubating::Float64
                            0.0, #time_recovering::Float64
                            0.0, #time_infectious_in_community::Float64
                            0.0, #test_sensitivity::Float64
                            0.0, #b1::Float64
                            0.0, #b2::Float64
                            0.0, #b3::Float64
                            0.0, #breakpoint::Float64
                            0.0, #beta_t::Float64
                            0.0, #beta_max::Float64
                            0.0, #beta_min::Float64
                            0.0, #k_inc::Float64
                            0.0, #k_r::Float64
                            0.0, #incubation_period::Float64
                            0.0, #recovery_period::Float64
                            false, #test_positive::Bool
                            false, #awaiting_test_result::Bool
                            0.0, #time_to_test_result::Float64
                            false, #discharged::Bool
                            false, #detected::Bool
                            false, #detected_clinical::Bool
                            false, #detected_test::Bool
                            false, # already_detected::Bool
                            false, #recorded::Bool
                            0.0,#FoI_in_community::Float64
                            0.0, #time_discharged::Float64
                            [0.0], #time_until_test::Array{Float64}
                            false #index_case::Bool
                            )



end
