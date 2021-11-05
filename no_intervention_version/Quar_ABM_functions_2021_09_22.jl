
using Distributions

# function definitions

function attempt_move!(agent::agent_type, destination, g, environments, model::AgentBasedModel)::Bool
  o = agent.pos
  d = destination
  if has_edge(g, o, d)
    move_agent!(agent, d, model)

    return true
    #println("moved between ", environments_2[o],  " and ",  environments_2[d])
  else
    #println("can't move between ", environments_2[o],  " and ",  environments_2[d])
    return false
  end
end
# checked


# initialise agents in specified model (note the below is for a GraphSpace)
function initialise_agents!(n_travellers, n_workers, model::AgentBasedModel, environments)
  for i in 1:(n_workers + n_travellers)
    if i <= n_workers
      add_agent!(environments["Hotel"], model)
      model[i].worker = true
      #println(model[i])
    else
      add_agent!(environments["Arrival"], model)
      model[i].traveller = true
      #println(model[i])
    end
  end
end
# checked

function initialise_work_schedule!(model::AgentBasedModel, environments, roster)
  days_of_the_week =  1:7;
  for i in keys(model.agents)
    if model[i].worker
      # assign work schedule
      shift_type = sample(rng, roster[1], roster[2])
      days_worked = shuffle(rng, days_of_the_week)
      days_worked[days_worked.> shift_type] .= 0
      days_worked[days_worked .> 0] .= 1
      workdays = days_of_the_week[days_worked .== 1]
      model[i].workdays = workdays
    end
  end
end
# checked

# check to make sure hotel is not understaffed:
function check_work_schedule(model::AgentBasedModel, n_min)::Bool
 days_of_the_week =  1:7
 workers_each_day = zeros(length(days_of_the_week))
  for i in keys(model.agents)
    if model[i].worker
      workers_each_day[model[i].workdays] .+= 1
    end
  end
  #println(workers_each_day)
  if minimum(workers_each_day) < n_min
    return false
  else
    return true
  end
end
# checked

# assigns uniform household sizes
function assign_households!(model::AgentBasedModel, n_per_hh, households)
  count = 0
  hh_id = 1
  households[hh_id] = household_type(hh_id, [])
  for i in keys(model.agents)
     if model[i].traveller
       model[i].hh_id = hh_id
       count += 1
       push!(households[hh_id].members, model[i].id)
       #println("assigned agent  $(model[i].id) to household $hh_id ")
       if count >= n_per_hh
         hh_id += 1
         households[hh_id] = household_type(hh_id, [])
         count = 0
       end
     end
   end

   # clean up empty households extras may be created due to workers and travellers
   # in teh same list.

   for hh_i in keys(households)
     if length(households[hh_i].members) == 0
       delete!(households, hh_i)
     end
   end
end
# checked

function vaccinate_all!(model::AgentBasedModel, pVac_Workers, pVac_Travellers)
  for i in keys(model.agents)
    if model[i].worker
      if rand(rng) < pVac_Workers
        model[i].vaccinated = true
      end
    elseif model[i].traveller
      if rand(rng) < pVac_Travellers
        model[i].vaccinated = true
      end
    end
  end
end
# checked


function vaccinate_infected_arrivals!(model::AgentBasedModel)
  for i in keys(model.agents)
    if model[i].traveller
      if environments_2[model[i].pos] == "Arrival"
        if model[i].infected && !model[i].vaccinated
          model[i].vaccinated = true
        end
      end
    end
  end
end
# checked

# infect arrivals and initialise infectivity and test sensitivity treajecotires

function infect_arrival!(agent::agent_type, # agent type
                         latent_period,
                         inc_dist, # incubation period distribution
                         r_dist,  # recovery period distribution
                         b_dist, # distribution of infectivity (gamma)
                         c_lims, # range for test sensitivity breakpoint
                         b1_lims, # range for max test sensitivity parameter b1
                         b2_lims, # range for test sensitivity growth parameter b2
                         b3_lims, # range for test sensitivity decay parameter b3
                         p_asymp # probability of becoming symptomatic after incubation
                         )

  agent.infected = true
  agent.index_case = true
  if rand(rng) < p_asymp
    agent.asymptomatic = true
  end
  while agent.incubation_period < latent_period
    agent.incubation_period = rand(rng, inc_dist)
  end
  agent.recovery_period = rand(rng, rec_dist)
  agent.beta_max = rand(rng, b_dist)
  agent.beta_min = agent.beta_max / Vmax

  # note the kinc and krec functions are in the 'disease' file
  agent.k_inc = compute_kinc(agent.incubation_period)
  agent.k_r = compute_krec(agent.recovery_period)

  q = cdf(inc_dist, agent.incubation_period)
  # set test sensitivity breakpoint - it's a function of incubation period
  agent.breakpoint = agent.incubation_period - ( q * (c_lims[2] - c_lims[1]))

  #set test sensitivity growth rate - also a function of incubation period
  agent.b2 = b2_lims[1] + (1.0 - q) * (b2_lims[2] - b2_lims[1])
  # might as well initialise the other test sensitivity terms too
  agent.b1 = b1_lims[1] + rand(rng) * (b1_lims[2] - b1_lims[1])
  agent.b3 = b3_lims[1] + rand(rng) * (b3_lims[2] - b3_lims[1])

  # determine infection status (timing) upon arrival
  # must be pre-symptomatic if not asymptomatic.

  if agent.asymptomatic
    agent.time_infected = rand(rng) * (agent.incubation_period + agent.recovery_period)
  else
    agent.time_infected = rand(rng) * agent.incubation_period
  end

  # determine current state based on time infected
  if agent.time_infected < latent_period
    agent.latent = true
    agent.incubating = true
    agent.time_incubating = agent.time_infected
  elseif agent.time_infected < agent.incubation_period
        agent.incubating = true
        agent.time_incubating = agent.time_infected
  else
         agent.recovering = true
         agent.time_incubating = agent.incubation_period
         agent.time_recovering = agent.time_infected - agent.incubation_period
  end
end
# checked

function infect_arrivals!(model::AgentBasedModel, pI_0, pI_Vac)
  n_arrivals_infected = 0
  for i in keys(model.agents)
    if model[i].traveller
      if environments_2[ model[i].pos ] == "Arrival"
        if (model[i].vaccinated && rand(rng) < pI_Vac) ||
           (!model[i].vaccinated && rand(rng) < pI_0)

           n_arrivals_infected += 1
            infect_arrival!(model[i], # agent type
                            latent_period,
                            inc_dist, # incubation period distribution
                            rec_dist,  # recovery period distribution
                            b_dist, # distribution of infectivity (gamma)
                            c_lims, # range for test sensitivity breakpoint
                            b1_lims, # range for max test sensitivity parameter b1
                            b2_lims, # range for test sensitivity growth parameter b2
                            b3_lims, # range for test sensitivity decay parameter b3
                            p_asymp # probability of becoming symptomatic after incubation)
                            )
            #println("arrival infected")
        end
      end
    end
  end

 println("infected $n_arrivals_infected arrivals")
end
# checked


function infect_agent!(agent::agent_type)

  agent.infected = true
  agent.time_infected = 0.0
  agent.latent = true
  agent.incubating = true
  agent.recovering = false
  agent.recovered = false
  agent.expressing_symptoms = false
  agent.time_infected = 0.0

  if rand(rng) < p_asymp
    agent.asymptomatic = true
  end
  while agent.incubation_period < latent_period
    agent.incubation_period = rand(rng, inc_dist)
  end
  agent.recovery_period = rand(rng, rec_dist)
  agent.beta_max = rand(rng, b_dist)
  agent.beta_min = agent.beta_max / Vmax

  # note the kinc and krec functions are in the 'disease' file
  agent.k_inc = compute_kinc(agent.incubation_period)
  agent.k_r = compute_krec(agent.recovery_period)

  q = cdf(inc_dist, agent.incubation_period)
  # set test sensitivity breakpoint - it's a function of incubation period
  #agent.breakpoint = c_lims[1] + q * (c_lims[2] - c_lims[1])
  agent.breakpoint = agent.incubation_period - ( q * (c_lims[2] - c_lims[1]))


  #set test sensitivity growth rate - also a function of incubation period
  agent.b2 = b2_lims[1] + (1.0 - q) * (b2_lims[2] - b2_lims[1])
  # might as well initialise the other test sensitivity terms too
  agent.b1 = b1_lims[1] + rand(rng) * (b1_lims[2] - b1_lims[1])
  agent.b3 = b3_lims[1] + rand(rng) * (b3_lims[2] - b3_lims[1])
end
# checked

function assign_workforce!(model::AgentBasedModel, day_of_week)::Float64
    n_present = 0.0
    for i in keys(model.agents)
          if model[i].worker
                if day_of_week in model[i].workdays
                #  println("working")
                  model[i].present_in_workplace = true
                  n_present += 1.0
                    if environments_2[model[i].pos] == "Community"
                      attempt_move!(model[i], environments["Hotel"], g, environments, model)
                    end
                else
                #  println("not working")
                  model[i].present_in_workplace = false
                    if environments_2[model[i].pos] == "Hotel"
                      attempt_move!(model[i], environments["Community"], g, environments, model)
                    end
                end
          end
    end

    return n_present
end
# checked

function compute_test_sensitivity!(agent::agent_type)
  #println("computing test sensitivity")
    agent.test_sensitivity = 0.0
    if (agent.infected || agent.recovered)

          tau = agent.time_infected - agent.breakpoint
          #println("computing tau: $tau" )

          if tau < 0.0
            agent.test_sensitivity = 1.0 / (1.0 + exp(-1.0 * (agent.b1 + agent.b2 * tau )))
            #println(1.0 / (1.0 + exp(-1.0 * (agent.b1 + agent.b2 * tau ))))

          else
            agent.test_sensitivity = 1.0 / (1.0 + exp(-1.0 * (agent.b1 + agent.b2 * tau + (-1.0 * agent.b2*agent.b3*tau) )))
            #println(1.0 / (1.0 + exp(-1.0 * (agent.b1 + agent.b2 * tau + agent.b2*agent.b3*tau ))))
          end

          #HACK: forcing test sensitivity of 0.5
          #agent.test_sensitivity = 0.5

    end
end
# checked

function check_test_sensitivity!(model::AgentBasedModel, t_sens, i_check)
  for i in keys(model.agents)
    if model[i].test_sensitivity > 0 && !(i in i_check)
      push!(t_sens, model[i].test_sensitivity)
      push!(i_check, i)
    end
  end
  if length(t_sens) > 0
   println("average test sensitivity: $(mean(t_sens))")
   println("tests performed: $(length(i_check))")
  end
end
# checked

function test_agent!(agent::agent_type, delay)

    compute_test_sensitivity!(agent)

    #println("test sensitivity: $(agent.test_sensitivity)")
      # test sensitivity can only be above 0 if agent is infected or recovered
        if rand(rng) < agent.test_sensitivity
          #println("test is positive")
          agent.test_positive = true
        else
          #println("test is negative")
          agent.test_positive = false
        end

    agent.time_to_test_result = delay
    agent.awaiting_test_result = true;
end
# checked

function test_workforce!(model::AgentBasedModel, delay)
    for i in keys(model.agents)
        if model[i].worker
            if model[i].present_in_workplace
               test_agent!(model[i], delay)
            end
        end
    end
end
# checked

function test_travellers!(model::AgentBasedModel, delay)
  for i in keys(model.agents)
      if model[i].traveller
        #println("time until traveller test: $(minimum(model[i].time_until_test)) ")
          if minimum(model[i].time_until_test) <= dt
             test_agent!(model[i], delay)
             #println("$(model[i].time_until_test)")
             #println("$(minimum(model[i].time_until_test))")
             #println("time until traveller $i test: $((model[i].time_until_test)) ")
             model[i].time_until_test[argmin(model[i].time_until_test)] = Inf
             #println("$(model[i].time_until_test)")
             #println("time until traveller $i test: $((model[i].time_until_test)) ")
          end
      end
  end
end
# checked

function assess_symptoms!(model::AgentBasedModel, iso_duration, iso_symptom_extension_flag)
    for i in keys(model.agents)
        if model[i].expressing_symptoms &&
          environments_2[model[i].pos] in ["Hotel",
                                           "Hotel Extension",
                                           "Home",
                                           "Home Extension",
                                           "MediHotel Extension",
                                           "Isolation"]

            if (!model[i].detected)
               model[i].detected = true
               if model[i].traveller
                 if model[i].discharged
                    model[i].discharged = false
                 end
               end

            elseif (environments_2[model[i].pos] == "Isolation" &&
                    !model[i].isolation_extended && iso_symptom_extension_flag)
               model[i].time_remaining_in_isolation = iso_duration
               model[i].isolation_extended = true
            end
        end
    end
end
# checked

function evaluate_workforce_tests!(model::AgentBasedModel)
  for i in keys(model.agents)
      if model[i].worker
        if (model[i].awaiting_test_result &&
           model[i].time_to_test_result <= dt)

           model[i].awaiting_test_result = false

           if model[i].test_positive
            model[i].detected = true
           end

        end
      end
  end
end
# checked

function evaluate_traveller_tests!(model::AgentBasedModel)
  for i in keys(model.agents)
      if model[i].traveller
          if (model[i].awaiting_test_result &&
              model[i].time_to_test_result <= dt)

             model[i].awaiting_test_result = false

             if model[i].test_positive
                model[i].detected = true
                if model[i].discharged
                   model[i].discharged = false
                end
             end

          end
      end
  end
end
# checked

function reinitialise_worker!(agent::agent_type, pVac_W, p_asymptomatic)
  # work schedule stays the same, id and position stay the same,
  # other parameters are reset:
  agent.worker = true
  agent.traveller = false
  agent.vaccinated = (rand(rng) < pVac_W)
  agent.infected = false
  agent.latent = false
  agent.incubating = false
  agent.recovering = false
  agent.recovered = false
  agent.asymptomatic = (rand(rng) < p_asymptomatic)
  agent.expressing_symptoms = false
  agent.time_infected = 0.0
  agent.time_incubating = 0.0
  agent.time_recovering = 0.0
  agent.time_infectious_in_community = 0.0
  agent.test_sensitivity = 0.0
  agent.b1 = 0.0
  agent.b2 = 0.0
  agent.b3 = 0.0
  agent.breakpoint = 0.0
  agent.beta_t = 0.0
  agent.beta_max = 0.0
  agent.beta_min = 0.0
  agent.k_inc = 0.0
  agent.k_r = 0.0
  agent.incubation_period = 0.0
  agent.recovery_period = 0.0
  agent.test_positive = false
  agent.time_to_test_result = 0.0
  agent.discharged = false
  agent.detected = false
  agent.recorded = false
  agent.FoI_in_community = 0.0
  agent.time_discharged = 0.0
  agent.time_until_test = [0.0]
  agent.awaiting_test_result = false;
end
# checked

function record_worker(agent::agent_type)
  linelist_W_entry[:exposure_days] = agent.time_infectious_in_community
  linelist_W_entry[:t_latent] = latent_period
  linelist_W_entry[:t_incubation] = agent.incubation_period
  linelist_W_entry[:t_post_incubation] = agent.recovery_period
  linelist_W_entry[:tested_positive] = convert(Int64, agent.test_positive)
  linelist_W_entry[:expressed_symptoms] = convert(Int64, (!agent.asymptomatic & (agent.time_recovering > 0)))
  linelist_W_entry[:time_discharged] = agent.time_discharged
  linelist_W_entry[:symptomatic] = convert(Int64, !agent.asymptomatic)
  linelist_W_entry[:vaccinated] = convert(Int64, agent.vaccinated)
  linelist_W_entry[:FoI_max] = agent.beta_max
  linelist_W_entry[:FoI_community] = agent.FoI_in_community
  agent.recorded = true
  push!(linelist_W, linelist_W_entry)
end
# checked

function discharge_workers!(model::AgentBasedModel, t, pVac_W, p_asymptomatic)
  for i in keys(model.agents)
      if model[i].worker
          if model[i].detected
             model[i].discharged = true
             model[i].time_discharged = t
             if model[i].recovered
               println("recovered worker recorded")
             end
             record_worker(model[i])
             reinitialise_worker!(model[i], pVac_W, p_asymptomatic)
          end
      end
  end
end
# checked

function compute_FoI!(agent::agent_type)
    if !agent.infected
        agent.beta_t = 0.0
      elseif agent.incubating && !agent.latent
          t_inc = agent.time_infected
          agent.beta_t = agent.beta_min * ( exp( agent.k_inc * t_inc) / Vmax)
        else
          t_r = agent.time_infected - agent.incubation_period
          agent.beta_t = agent.beta_max * exp(agent.k_r * t_r)
    end
    #plateau implemented as a cutoff.
    if agent.beta_t > agent.beta_max
      agent.beta_t = agent.beta_max
    end
end
# checked

function step_infection!(agent::agent_type, time_increment)::Float64
    exposure_days = 0.0
    #update infection step

    infectious_period = agent.incubation_period + agent.recovery_period

    if agent.time_infected > infectious_period

           agent.infected = false
           agent.latent = false
           agent.incubating = false
           agent.recovering = false
           agent.recovered = true
           agent.expressing_symptoms = false

           if agent.worker
             println("worker recovered!")
             agent.detected = true # a hack to make sure
             # recovered workers are swapped out with susceptible ones
             # and recorded in the line list.
           end



    elseif (agent.time_infected > latent_period) &&
           (agent.time_infected < agent.incubation_period)

           agent.latent = false
           agent.incubating = true
           agent.recovering = false
           agent.recovered = false
           agent.expressing_symptoms = false

    elseif agent.time_infected > agent.incubation_period

           agent.latent = false
           agent.incubating = false
           agent.recovering = true
           agent.recovered = false

           if !agent.asymptomatic
              agent.expressing_symptoms = true
           end

    end

    # iterate clocks related to infection
    if agent.infected
      agent.time_infected += time_increment
      if ( (environments_2[agent.pos] == "Community" || agent.worker ) &&  !agent.latent) ||
         ( agent.traveller && !agent.compliant && !agent.latent )
        agent.FoI_in_community += (agent.beta_t * time_increment)
        agent.time_infectious_in_community += time_increment
        exposure_days += time_increment
      end
    end

    if agent.incubating
       agent.time_incubating += time_increment
    elseif agent.recovering
      agent.time_recovering += time_increment
    end

    # compute force of infection for next iteration
    compute_FoI!(agent)

    return exposure_days
end
# checked

function step_workers!(model::AgentBasedModel, time_increment::Float64)::Tuple{Float64, Float64, Float64}

    exposure_days = 0.0
    exposure_days_uv = 0.0
    exposure_days_vac = 0.0

    for i in keys(model.agents)
        if model[i].worker
            if model[i].infected
              exposure_days_i = step_infection!(model[i], time_increment)
              exposure_days += exposure_days_i
              if model[i].vaccinated
                exposure_days_vac += exposure_days_i
              else
                exposure_days_uv += exposure_days_i
              end
            end
          # iterate clocks related to testing
            if model[i].awaiting_test_result
               model[i].time_to_test_result -= dt
            end
        end
    end

    return exposure_days, exposure_days_uv, exposure_days_vac
end
# checked

function isolate_travellers!(model::AgentBasedModel,
                             iso_duration,
                             test_schedule_ext,
                             medi_hotel_flag,
                             home_ext_compliance_rate)::Tuple{Float64, Float64}
  PCR_detections = 0.0
  clinical_detections = 0.0

  for i in keys(model.agents)
      if model[i].traveller
          if (model[i].detected && environments_2[model[i].pos] in ["Home",
                                                                    "Hotel",
                                                                    "MediHotel Extension",
                                                                    "Hotel Extension",
                                                                    "Home Extension"] )
              attempt_move!(model[i], environments["Isolation"],g,environments, model)
              if model[i].discharged
                println("trying to isolate discharged traveller $i after $(model[i].time_in_quar) days in quar")
                println("agent expressing symptoms? $(model[i].expressing_symptoms)")
              end
              # NOTE: if the person expressed symptoms right after returning the positive
              # test, it's counted as a clinical detection.
              if model[i].expressing_symptoms
                clinical_detections += 1.0
              else
                 PCR_detections += 1.0
              end


              model[i].detected = false
              model[i].time_remaining_in_isolation = iso_duration
              model[i].time_until_test = ones(size(model[i].time_until_test)) .* Inf
              #println(["isolated agent $i, now in" environments_2[model[i].pos]])
              # extend quarantine for contacts
              # NOTE: if a contact is already detected, they will get moved when
              # their turn comes - if they've already been moved to iso, the
              # attempt_move function will return false.
              hh = model[i].hh_id
              hh_members = households[hh].members

              for j in hh_members

                  if !model[j].detected

                      moved = false
                      if model[j].discharged && !model[j].quarantine_extended
                        model[j].discharged = false
                      end

                      if medi_hotel_flag
                        previous = environments_2[model[j].pos]
                         moved = attempt_move!(model[j], environments["MediHotel Extension"],g,environments, model)
                           if moved
                             model[j].p_compliance = 1.0
                             if model[j].discharged
                               println("trying to move discharged traveller $j from $(previous) to MediHotel after $(model[j].time_in_quar) days in quar")
                             end
                           end
                        elseif environments_2[model[j].pos] == "Home"
                          moved = attempt_move!(model[j], environments["Home Extension"],g,environments, model)
                            if moved
                            model[j].p_compliance = home_ext_compliance_rate
                            end
                        elseif environments_2[model[j].pos] == "Hotel"
                          moved = attempt_move!(model[j], environments["Hotel Extension"],g,environments, model)
                            if moved
                              model[j].p_compliance = 1.0
                            end
                      end

                      if moved
                        model[j].time_until_test = copy(test_schedule_ext)
                        model[j].quarantine_extended = true
                      #else
                        #println("tried to move from $(environments_2[model[j].pos]) to extension")
                      end
                  end
              end
          end
      end
  end

  return PCR_detections, clinical_detections
end
# checked

function discharge_travellers!(model::AgentBasedModel, t, n_new_hh, n_agents_to_add)::Tuple{Float64, Float64}

    n_discharged_t = 0.0
    n_discharged_infected_t = 0.0

    hh_of_discharged_agents = []

    for i in keys(model.agents)
        if model[i].traveller
            if model[i].discharged && model[i].time_discharged == 0.0
               model[i].time_discharged = t

               tag = environments_2[model[i].pos]

               if tag == "Community"
                 println("trying to discharge twice")
               end

               flag = attempt_move!(model[i], environments["Community"], g, environments, model)
               n_discharged_t += 1.0
               if model[i].infected
                 n_discharged_infected_t += 1.0
               end
               #=if tag == "Isolation"
                  if model[i].infected
                    println(["isolated, infected traveller $i moved to community:  $flag"])
                  elseif model[i].recovered
                    println(["isolated, recovered traveller $i moved to community:  $flag"])
                  else
                    println(["isolated, susceptible traveller $i moved to community:  $flag"])
                  end
               end=#

               #note the below may need to be reinstated.
               #model[i].discharged = false
               if !(model[i].hh_id in hh_of_discharged_agents)
                  push!(hh_of_discharged_agents, model[i].hh_id)
               end

             end
         end
     end


     for hh_i in hh_of_discharged_agents
      # iterate through household members and see if they're all in the community
         n_hh = length(households[hh_i].members)
         n_in_community = 0

         for j in households[hh_i].members
           if environments_2[model[j].pos] == "Community"
             n_in_community += 1
           end
         end

         if n_in_community == n_hh
           n_new_hh .+= 1
           n_agents_to_add .+= n_in_community
         end
     end

 return n_discharged_t, n_discharged_infected_t
end
# checked


function initialise_arrivals!(model::AgentBasedModel,
                              pVac_T,
                              pI_Vac,
                              pI_0,
                              n_new_hh,
                              n_new_agents,
                              test_schedule)::Tuple{Float64, Float64, Float64}
    n_arrivals_infected = 0.0
    n_arrivals_infected_uv = 0.0
    n_arrivals_infected_vac = 0.0
    # generate the new agents:
    new_ids = []
    for i in 1:n_new_agents
      new_id = nextid(model)
      push!(new_ids, new_id)
      add_agent!(environments["Arrival"], model)
      model[new_id].traveller = true
    end

    #assign households to new agents:

    n_per_hh = n_new_agents / n_new_hh # should be constant in this version

    if length(households) == 0
      hh_id = 1
        else
        hh_id = maximum(keys(households)) + 1
    end
    count = 0
    households[hh_id] = household_type(hh_id, [])
    for i in new_ids

        if count >= n_per_hh
           hh_id += 1
           households[hh_id] = household_type(hh_id, [])
           count = 0
        end

         model[i].hh_id = hh_id
         count += 1
         push!(households[hh_id].members, model[i].id)
         #println("assigned agent  $(model[i].id) to household $hh_id ")
    end

    # vaccinate, infect, and move new arrivals
    for i in new_ids

        if rand(rng) < pVac_T
           model[i].vaccinated = true
        end

        if (model[i].vaccinated && rand(rng) < pI_Vac) ||
           (!model[i].vaccinated && rand(rng) < pI_0)

            infect_arrival!(model[i], # agent type
                            latent_period,
                            inc_dist, # incubation period distribution
                            rec_dist,  # recovery period distribution
                            b_dist, # distribution of infectivity (gamma)
                            c_lims, # range for test sensitivity breakpoint
                            b1_lims, # range for max test sensitivity parameter b1
                            b2_lims, # range for test sensitivity growth parameter b2
                            b3_lims, # range for test sensitivity decay parameter b3
                            p_asymp # probability of becoming symptomatic after incubation)
                            )
            n_arrivals_infected += 1.0
            if model[i].vaccinated
              n_arrivals_infected_vac += 1.0
            else
              n_arrivals_infected_uv += 1.0
            end
            #println("arrival infected")
        end

        #attempt_move!(model[i], environments["Hotel"], g, environments)
        #model[i].time_until_test = copy(test_schedule)

    end

  #=n_total = length(new_ids)
  if n_total > 0
  println("infected $n_arrivals_infected out of $n_total arrivals ")
  end=#

 return n_arrivals_infected, n_arrivals_infected_uv, n_arrivals_infected_vac
end
# checked

function move_arrivals_to_home_quarantine!(model::AgentBasedModel, test_schedule)
    for i in keys(model.agents)
        if environments_2[model[i].pos] == "Arrival"
          attempt_move!(model[i], environments["Home"], g, environments, model)
          model[i].time_until_test = copy(test_schedule)
        end
    end
end
# checked

function move_arrivals_to_hotel_quarantine!(model::AgentBasedModel, test_schedule)
    for i in keys(model.agents)
        if environments_2[model[i].pos] == "Arrival"
          attempt_move!(model[i], environments["Hotel"], g, environments, model)
          model[i].time_until_test = copy(test_schedule)
          #println("FUCK")
        end
    end
end
# checked


function record_and_delete!(model::AgentBasedModel)
   # determine who to kill
   # iterate through households and evaluate removal criteria
   agents_to_kill = []
   households_to_kill = []
   for hh_i in keys(households)
       n_hh = length(households[hh_i].members)
       n_kill = 0
       for i in households[hh_i].members
           if (environments_2[model[i].pos] == "Community" &&
               model[i].infected == false)
               n_kill += 1
           end
       end
       if n_kill == n_hh
           push!(households_to_kill, hh_i )
           agents_to_kill = vcat(agents_to_kill, households[hh_i].members)
       end
   end

   for i in agents_to_kill
      if model[i].infected || model[i].recovered
         record_traveller(model[i])
      end

       kill_agent!(model[i], model)
   end

   for hh_i in households_to_kill
     delete!(households, hh_i)
   end
end
# checked

function record_traveller(agent::agent_type)
  linelist_T_entry[:exposure_days] = agent.time_infectious_in_community
  linelist_T_entry[:days_in_quar] = agent.time_in_quar
  linelist_T_entry[:days_in_extended_quar] = agent.time_in_extended_quar
  linelist_T_entry[:days_in_iso] = agent.time_in_isolation
  linelist_T_entry[:t_latent] = latent_period
  linelist_T_entry[:t_incubation] = agent.incubation_period
  linelist_T_entry[:t_post_incubation] = agent.recovery_period
  linelist_T_entry[:time_discharged] = agent.time_discharged
  linelist_T_entry[:index_case] = convert(Int64, agent.index_case)
  linelist_T_entry[:symptomatic] = convert(Int64, !agent.asymptomatic)
  linelist_T_entry[:vaccinated] = convert(Int64, agent.vaccinated)
  linelist_T_entry[:compliance] = agent.p_compliance
  linelist_T_entry[:FoI_max] = agent.beta_max
  linelist_T_entry[:FoI_community] = agent.FoI_in_community
  agent.recorded = true
  push!(linelist_T, linelist_T_entry)
end
# checked

function step_travellers!(model::AgentBasedModel,
                          quar_duration,
                          quar_extension,
                          iso_duration,
                          time_increment)::Tuple{Float64, Float64, Float64}
    exposure_days = 0.0
    exposure_days_uv = 0.0
    exposure_days_vac = 0.0

    for i in keys(model.agents)
        if model[i].traveller
            if model[i].infected
                exposure_days_i = step_infection!(model[i], time_increment)
                exposure_days += exposure_days_i
                if model[i].vaccinated
                  exposure_days_vac += exposure_days_i
                else
                  exposure_days_uv += exposure_days_i
                end
            end
            # iterate clocks related to testing
            if model[i].awaiting_test_result
               model[i].time_to_test_result -= time_increment
            end

            for k in 1:length(model[i].time_until_test)
              model[i].time_until_test[k] -= time_increment
            end

            # remember to update 'discharged' status
            #if environments_2[model[i].pos] in ["Hotel", "Home"]
            if (environments_2[model[i].pos] == "Hotel" || environments_2[model[i].pos] == "Home" )
                model[i].time_in_quar += time_increment
                if model[i].time_in_quar > (quar_duration - time_increment) && !model[i].detected
                   model[i].discharged = true
                end

            elseif environments_2[model[i].pos] in ["Hotel Extension",
                                                    "MediHotel Extension",
                                                    "Home Extension"]
                  model[i].time_in_extended_quar += time_increment
                  if model[i].time_in_extended_quar > (quar_extension - time_increment) && !model[i].detected
                    model[i].discharged = true
                  end

            elseif environments_2[model[i].pos] == "Isolation"

                    model[i].time_in_isolation += time_increment
                    model[i].time_remaining_in_isolation -= time_increment
                    #println("$i has been in iso for: $(model[i].time_in_isolation)")
                    #println("$i has : $(model[i].time_remaining_in_isolation) remaining")

                    #if model[i].time_remaining_in_isolation < -1.0
                    #=if ((model[i].time_in_isolation) + (model[i].time_remaining_in_isolation)) > 11
                      println("$i total iso time of:
                      $(model[i].time_in_isolation) + $(model[i].time_remaining_in_isolation) = $((model[i].time_in_isolation) + (model[i].time_remaining_in_isolation))")
                    end=#



                    if model[i].time_remaining_in_isolation <= dt
                      model[i].discharged = true
                      #println("$i has : $(model[i].time_remaining_in_isolation) remaining")
                      #println("$i should be discharged.")
                      #println("discharging from isolation: $(model[i].infected)" )
                    end

            end
          end
    end

    return exposure_days, exposure_days_uv, exposure_days_vac
end
# checked


function compute_transmission!(model::AgentBasedModel,
                               VEi,
                               VEs,
                               FoI_fac_Hotel_TT_same,
                               FoI_fac_Hotel_TT_dif,
                               FoI_fac_Hotel_TW,
                               FoI_fac_Hotel_WW,
                               n_workers_present,
                               time_increment)::Float64
 n_transmissions = 0.0
 n_travellers_in_hotel =
 convert(Float64, length(agents_in_position(environments["Hotel"], model))) +
 convert(Float64, length(agents_in_position(environments["Hotel Extension"], model)))
 - n_workers_present

    for i in keys(model.agents)

        if model[i].infected &&
          !model[i].latent &&
          !(environments_2[model[i].pos] in ["Arrival", "Community", "Isolation", "MediHotel Extension"])

            for j in keys(model.agents)
                #if model[i].pos == model[j].pos
                if !model[j].infected && !model[j].recovered

                    env_fac = 0.0
                    pop_fac = 1.0
                    vac_fac = (1.0 - (VEs * convert(Float64, model[i].vaccinated) )  ) *
                              (1.0 - (VEi * convert(Float64, model[j].vaccinated) )  )

                    # check if agents are located in interacting environments
                  #  if (environments_2[model[i].pos] in ["Hotel", "Hotel Extension"] ) &&
                  #     (environments_2[model[j].pos] in ["Hotel", "Hotel Extension"] )


                  if (environments_2[model[i].pos] == "Hotel" || environments_2[model[i].pos] == "Hotel Extension" ) &&
                     (environments_2[model[j].pos] == "Hotel" || environments_2[model[j].pos] == "Hotel Extension" )
                        # transmit to household members
                        if (model[i].traveller && model[j].traveller && model[i].hh_id == model[j].hh_id)

                           n_HH = convert(Float64, length(households[model[i].hh_id].members)) - 1.0
                           pop_fac = 1.0 / n_HH
                           env_fac = FoI_fac_Hotel_TT_same

                                  # transmit between all other agents in hotel
                                  # travellers and workers:

                        elseif (model[i].traveller && model[j].traveller && model[i].hh_id != model[j].hh_id)

                            n_HH = length(households[model[i].hh_id].members)
                            n_hotel = n_travellers_in_hotel - n_HH
                            pop_fac = 1.0 / n_hotel
                            env_fac = FoI_fac_Hotel_TT_dif

                        elseif (model[i].worker && model[j].traveller)

                            n_hotel = n_travellers_in_hotel
                            pop_fac = 1.0 / n_hotel
                            env_fac = FoI_fac_Hotel_TW

                        elseif (model[i].traveller &&
                               (model[j].worker && model[j].present_in_workplace ) )

                            n_hotel = n_workers_present
                            pop_fac = 1.0 / n_hotel
                            env_fac = FoI_fac_Hotel_TW

                        elseif (model[i].worker &&
                               (model[j].worker && model[j].present_in_workplace ) )

                            n_hotel = n_workers_present - 1.0
                            pop_fac = 1.0 / n_hotel
                            env_fac = FoI_fac_Hotel_WW
                        end

                    elseif environments_2[model[i].pos] in ["Home", "Home Extension"] &&
                             environments_2[model[j].pos] in ["Home", "Home Extension"]
                        # transmit between household members only
                        if (model[i].traveller &&
                            model[i].hh_id == model[j].hh_id)

                                  n_HH = convert(Float64, length(households[model[i].hh_id].members)) - 1.0
                                  pop_fac = 1.0 / n_HH
                                  env_fac = FoI_fac_Hotel_TT_same
                        end
                    end

                    beta_ij = model[i].beta_t * time_increment * pop_fac * env_fac * vac_fac

                    p_infect = 1.0 - exp((-1.0 * beta_ij))

                  #=  if model[j].worker
                      println("b infect worker: $p_infect")
                    end=#

                    #=if model[j].traveller
                      println("p infect traveller: $p_infect")
                    end=#


                    if rand(rng) < p_infect
                      infect_agent!(model[j])
                      n_transmissions += 1.0
                        #=if model[j].worker
                          println("worker infected")
                        end=#
                        #=if model[j].traveller
                          println("traveller infected")
                        end=#
                    end
                end
            end
        end
    end

    return n_transmissions
end
# checked

function evaluate_compliance!(model)
  for i in keys(model.agents)
    if model[i].traveller &&
      environments_2[model[i].pos] in ["Home", "Home Extension"]

       model[i].compliant = (rand(rng) < model[i].p_compliance)
       #println("agent $i compliant? : $(model[i].compliant)")
       #println("compliance probability: $(model[i].p_compliance)")
     else
       model[i].compliant = true
     end
   end
end
# checked

function assign_compliance_rate_to_arrivals!(model::AgentBasedModel, compliance_rate)
  for i in keys(model.agents)
    #model[i].p_compliance = 1.0
    if model[i].traveller &&
       environments_2[model[i].pos] == "Arrival"
       model[i].p_compliance = compliance_rate
    end
   end
end
# checked

function check_agents!(model::AgentBasedModel, t)

  for i in keys(model.agents)

    if model[i].time_in_isolation > 0.1 &&
      environments_2[model[i].pos] in ["Home", "Hotel", "MediHotel Extension", "Home Extension", "Hotel Extension"]
      println("***** $i bad move ******")
    end

    if model[i].time_in_isolation > 0.1 &&
      !(environments_2[model[i].pos] in ["Isolation", "Community"])
      println("***** $i bad move ******")
    end

    if model[i].p_compliance == 1.0 &&
      #(environments_2[model[i].pos] != "MediHotel Extension") &&
      #(environments_2[model[i].pos] != "Community") &&
      #(environments_2[model[i].pos] != "Isolation") &&
      !model[i].quarantine_extended &&
      model[i].traveller#!(model[i].quarantine_extended)
      println("****compliance changed weirdly for agent $i at time $t *****")
      println("****$(model[i].pos) => $(environments_2[model[i].pos])*****")
    end


  end
end
# checked



            # if the above conditions are satisfied, transmission is possible under some circumstances
