abstract type AbstractMyResult end 

mutable struct MyODEResult <: AbstractMyResult
    points::Vector{Float64}
    values::Vector{Float64}
    a_tol::Vector{Float64}
    r_tol::Vector{Float64}

    function MyODEResult(x,y)
        a_tol=[NaN for _ in eachindex(x)]
        r_tol=[NaN for _ in eachindex(x)]
        new{}(x, y, a_tol, r_tol)
    end

    function MyODEResult(x, y, a_tol)
        r_tol=[NaN for _ in eachindex(x)]
        new{}(x, y, a_tol, r_tol)
    end
    
    function MyODEResult(x,y, a_tol, r_tol)
        new{}(x, y, a_tol, r_tol)
    end
end
