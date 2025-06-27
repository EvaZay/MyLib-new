abstract type AbstractMyProblem end 

mutable struct MyODEProblem{Method, FType, RType} <: AbstractMyProblem
    funcs::FType
    x::Float64
    # a_tol::Vector{Float64}
    # r_tol::Vector{Float64}
    results::Vector{Float64}
    st::Float64
    coeffs1::Vector{Float64}
    coeffs2::Vector{Float64}
    coeffs3::Vector{Float64}
    coeffs4::Vector{Float64}
    coeffs5::Vector{Float64}
    cash::Vector{Float64}
    function MyODEProblem(funcs::Vector{<:Function}, method::String, x0::Real, initial_conditions::Vector{<:Real}, st::Real)
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

        N=length(funcs)

        coeffs1=zeros(Float64, N)
        coeffs2=zeros(Float64, N)
        coeffs3=zeros(Float64, N)
        coeffs4=zeros(Float64, N)
        coeffs5=zeros(Float64, N)

        # a_tol=zeros(Float64, n)

        # r_tol=zeros(Float64, n)
        
        new{Method, typeof(funcs), typeof(initial_conditions)}(funcs, x0,  initial_conditions, st, coeffs1, coeffs2, coeffs3, coeffs4, coeffs5)
    end

    function MyODEProblem(funcs::FType, method::String, x0::Real, initial_conditions::Real, st::Real) where FType<:Function
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
        
        coeffs1=0
        coeffs2=0
        coeffs3=0
        coeffs4=0
        coeffs5=0

        # a_tol=zeros(Float64, n)

        # r_tol=zeros(Float64, n)
        
        new{Method, typeof([funcs]), typeof(initial_conditions)}([funcs], x0, [initial_conditions], st, [coeffs1],  [coeffs2], [coeffs3], [coeffs4], [coeffs5])
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