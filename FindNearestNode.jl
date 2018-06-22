# Function to find the mesh node that is closest to the requested location

function FindNearestNode(xin, yin, X, Y)
    nseis = length(xin)
    dist = zeros(nseis, 1)
    iglob = zeros(Int32, nseis, 1)

    for k = 1:nseis
        dist[k], iglob[k] = findmin( (X - xin[k]).^2 + (Y - yin[k]).^2 )
    end

    dist = sqrt.(dist)
    xout = X[iglob]
    yout = Y[iglob]

    return xout, yout, iglob, dist

end
