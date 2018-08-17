# trying out parallel scripts

@everywhere function count_heads(n)
    c::Int = 0

    for i = 1:n

        c += rand(1)[1]
    end
    c
end
