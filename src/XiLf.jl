# Calculate XiLf used in computing the timestep

function XiLfFunc(s::space_parameters, t::time_parameters, muMax)

    hcell = s.LX/(s.FltNglob-1)
    Ximax = 0.5
    Xithf = 1

    for j = 1:s.FltNglob

        # Compute time restricting parameters
        expr1 = -()
