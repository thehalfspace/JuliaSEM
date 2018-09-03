#################################
# READ OUTPUT FROM SIMULATION
#################################

using Serialize
open("data.out") do f
    global Op, sim_time
    Op = deserialize(f)
    sim_time = deserialize(f)
end
