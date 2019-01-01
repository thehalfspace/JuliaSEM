using Distributed
using LinearAlgebra

addprocs(4)

@everywhere function compute(i)
    pp = i*i*i*i/sin(log(i*i*i*i*i))
    return pp
end


function dopar()
    nt = 2000
    for it = 1:nt
        @sync for pid in workers()
            @async remotecall_fetch(compute,pid,it)
        end
    end
end


function dopar2()
    nt = 2000
    @sync @distributed for it = 1:nt
        compute(it)
    end
end
