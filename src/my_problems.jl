abstract type AbstractMyProblem end 

mutable struct MyODEProblem{Method, FType, RType} <: AbstractMyProblem
    funcs::FType
    x::Vector{Float64}
    initial_conditions::Vector{Float64}
    a_tol::Vector{Float64}
    r_tol::Vector{Float64}
    results::RType
    coeffs1::Vector{Float64}
    coeffs2::Vector{Float64}
    coeffs3::Vector{Float64}
    coeffs4::Vector{Float64}
    coeffs5::Vector{Float64}
    function MyODEProblem(funcs::Vector{<:Function}, method::String, x0::Real, initial_conditions::Vector{<:Real}, n::Int64)
        if method == "Euler"
            Method = :Euler
        end
        if method == "Runge-Kutta 2"
            Method = :RungeKutta2
        end
        if method == "Runge-Kutta 3"
            Method = :RungeKutta3
        end
        if method == "Runge-Kutta 4"
            Method = :RungeKutta4
        end
        if method == "Runge-Kutta 5"
            Method = :RungeKutta5
        end

        x=zeros(Float64, n+1)
        x[1]=x0

        N=length(funcs)
        results = [zeros(Float64, n+1) for _ in 1:N]
        for i=1:N
            results[i][1]=initial_conditions[i]
        end

        coeffs1=zeros(Float64, N)
        coeffs2=zeros(Float64, N)
        coeffs3=zeros(Float64, N)
        coeffs4=zeros(Float64, N)
        coeffs5=zeros(Float64, N)

        a_tol=zeros(Float64, n)

        r_tol=zeros(Float64, n)
        
        new{Method, typeof(funcs), typeof(results)}(funcs, x, initial_conditions, a_tol, r_tol, results, coeffs1, coeffs2, coeffs3, coeffs4, coeffs5)
    end

    function MyODEProblem(funcs::FType, method::String, x0::Real, initial_conditions::Real, n::Int64) where FType<:Function
        if method == "Euler"
            Method = :Euler
        end
        if method == "Runge-Kutta 2"
            Method = :RungeKutta2
        end
        if method == "Runge-Kutta 3"
            Method = :RungeKutta3
        end
        if method == "Runge-Kutta 4"
            Method = :RungeKutta4
        end
        if method == "Runge-Kutta 5"
            Method = :RungeKutta5
        end

        x=zeros(Float64, n+1)
        x[1]=x0

        
        results = [zeros(Float64, n+1)]
        
        results[1][1]=initial_conditions
        
        coeffs1=0
        coeffs2=0
        coeffs3=0
        coeffs4=0
        coeffs5=0

        a_tol=zeros(Float64, n)

        r_tol=zeros(Float64, n)
        
        new{Method, typeof([funcs]), typeof(results)}([funcs], x, [initial_conditions], a_tol, r_tol, results, [coeffs1],  [coeffs2], [coeffs3], [coeffs4], [coeffs5])
    end
end

mutable struct MyIntegretionProblem{Method} <: AbstractMyProblem
    method::String
    accum::Float64
    st::Float64

    function MyIntegretionProblem(method::String, ic::Real, st::Real) 
        if method == "Forward Euler"
            Method = :ForwardEuler
        end
        if method == "Backward Euler"
            Method = :BackwardEuler
        end
        if method == "Trapezoidal"
            Method = :Trapezoidal
        end
        if method == "Simpsons"
            Method = :Simpsons
        end

        new{Method}(method, ic, st)
    end
end

struct MyInterpolationProblem <: AbstractMyProblem

end