#################################
# SOME ROUGH SCRIPTS
# FOR TRYING OUT STUFF
#################################

using StatsBase
using PyPlot

# Plot sliprates for each event with depth
function test1(tStart, time_, SlipVel)

    for i = 1:length(SlipVel[1,:])

        if 
