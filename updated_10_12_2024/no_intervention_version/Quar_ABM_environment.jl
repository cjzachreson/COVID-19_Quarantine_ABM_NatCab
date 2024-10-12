# define environment
using Graphs
using StatsBase

households = Dict()

environments =
Dict(
     [("Arrival", 1),
     ("Hotel", 2) ,
     ("Home", 3),
     ("MediHotel Extension", 4),
     ("Hotel Extension", 5),
     ("Home Extension", 6),
     ("Isolation", 7),
     ("Community", 8)])
# reverse mapping
environments_2 = Dict(value => key for (key, value) in environments)

# specify possible transitions (this defines the basic pathway structure)
edges = [[("Arrival", "Hotel"), 1],
         [("Arrival", "Home"), 1],
         [("Arrival", "Community"), 0],
         [("Hotel", "Home"), 0],
         [("Hotel", "MediHotel Extension"), 1],
         [("Hotel", "Hotel Extension"), 1],
         [("Hotel", "Isolation"), 1],
         [("Hotel", "Community"), 1],
         [("Community", "Hotel"), 1], # workers can do this, travellers can't TODO enforce this
         [("Home", "MediHotel Extension"), 1],
         [("Home", "Home Extension"), 1],
         [("Home", "Isolation"), 1],
         [("Home", "Community"), 1],
         [("MediHotel Extension", "Isolation"), 1],
         [("MediHotel Extension", "Community"), 1],
         [("Home Extension", "Isolation"), 1],
         [("Home Extension", "Community"), 1],
         [("Hotel Extension", "Isolation"), 1],
         [("Hotel Extension", "Community"), 1],
         [("Isolation", "Community"), 1]]

 # initialise network representation of
 #environment with the right number of nodes
 ne = length(environments);
 g = Graphs.SimpleDiGraph(ne)

  # add the edges specified in the list above:
 for e in edges
   #println(e[2])
   if e[2] == 1
   Graphs.add_edge!(g, environments[e[1][1]], environments[e[1][2]])
   end
 end

 # Space is a digraph representing the transitions between different
 #quarantine environemnts
 #(arrival, hotel, home, medi-hotel extension, isolation, community).
 Space = Agents.GraphSpace(g)

#specify rosters - these are the work schedules used by worker agents to decide
# which days they're at work and which days they're in the community.
# each roster specifies the number of days worked per week, and a propprtion
# of the workforce working that number of days per week
rosters = Dict([
               (1, [[7; 5; 3; 1],AnalyticWeights([1.0; 0.0; 0.0; 0.0])]),
               (2, [[7; 5; 3; 1],AnalyticWeights([0.0; 1.0; 0.0; 0.0])]),
               (3, [[7; 5; 3; 1],AnalyticWeights([0.0; 0.6; 0.4; 0.0])]),
               (4, [[7; 5; 3; 1],AnalyticWeights([0.0; 0.6; 0.3; 0.1])])
              ])
