abstract type AbstractMyProblem end 

"""
Задача численного решения дифференциальных уравнений

    Поля структуры:
    funcs -- функция или финкции, которую/-ые необходимо продефференцировать
    x -- переменная x
    results -- массив результатов на каждом шаге st дифференцирования
    st -- шаг дифференцирования
    coeffs1 -- первые коэффициенты методов для решения задачи
    coeffs2 -- вторые коэффициенты методов для решения задачи
    coeffs3 -- третьи коэффициенты методов для решения задачи
    coeffs4 -- четвертые коэффициенты методов для решения задачи
    coeffs5 -- пятые коэффициенты методов для решения задачи

    Конструктор для решения ОДУ
    MyODEProblem(funcs::FType, method::String, x0::Real, initial_conditions::Real, st::Real) where FType<:Function

        funcs -- функция, которую необходимо продефференцировать
        method -- метод, которым надо решить задачу
        x0 -- начальное условие для x
        initial_conditions -- начальное условие
        st -- шаг дифференцирования


    Конструктор для решения системы ОДУ
    MyODEProblem(funcs::Vector{<:Function}, method::String, x0::Real, initial_conditions::Vector{<:Real}, st::Real)

        funcs -- функции, которые необходимо продефференцировать
        method -- метод, которым надо решить задачу
        x0 -- начальное условие для x
        initial_conditions -- начальные условия
        st -- шаг дифференцирования

"""
const ODEMethodsForODE = Dict("Euler" => :Euler, "Runge-Kutta 2" => :RungeKutta2, 
    "Runge-Kutta 3" => :RungeKutta3, "Runge-Kutta 4" => :RungeKutta4)

mutable struct MyODEProblem{Method, FType, RType} <: AbstractMyProblem
    funcs::FType
    x::Float64
    a_tol::Vector{Float64}
    r_tol::Vector{Float64}
    results::Vector{Float64}
    st::Float64
    coeffs1::Vector{Float64}
    coeffs2::Vector{Float64}
    coeffs3::Vector{Float64}
    coeffs4::Vector{Float64}

    function MyODEProblem(funcs::Vector{<:Function}, method::String, x0::Real, 
        initial_conditions::Vector{<:Real}, st::Real, n::Int64)

        haskey(ODEMethodsForODE, method) || 
            throw(ArgumentError("Неподдерживаемый метод решения дифференциальных уравнений"))

        isempty(funcs) && throw()
        isempty(initial_conditions) && throw()
        length(funcs) != length(initial_conditions) && throw()
        st <= 0 && throw()
        n <= 0 && throw()

        N = length(funcs)
        coeffs1 = zeros(Float64, N)
        coeffs2 = zeros(Float64, N)
        coeffs3 = zeros(Float64, N)
        coeffs4 = zeros(Float64, N)

        a_tol = zeros(Float64, n+1)

        r_tol = zeros(Float64, n+1)
        
        new{ODEMethodsForODE[method], typeof(funcs), typeof(initial_conditions)}(
            funcs, x0, a_tol, r_tol, initial_conditions, st, coeffs1, coeffs2, coeffs3, coeffs4)
    end

    function MyODEProblem(funcs::FType, method::String, x0::Real, initial_conditions::Real, st::Real, n::Int64) where FType<:Function

        haskey(ODEMethodsForODE, method) || 
            throw(ArgumentError("Неподдерживаемый метод решения дифференциальных уравнений"))

        st <= 0 && throw()
        n <= 0 && throw()
        
        coeffs1 = 0
        coeffs2 = 0
        coeffs3 = 0
        coeffs4 = 0

        a_tol = zeros(Float64, n + 1)

        r_tol = zeros(Float64, n + 1)
        
        new{ODEMethodsForODE[method], typeof([funcs]), typeof(initial_conditions)}([funcs], x0, a_tol, r_tol, [initial_conditions], st, [coeffs1],  [coeffs2], [coeffs3], [coeffs4])
    end
end


"""
Задача численного интегрирования

    Поля структуры:
    method -- метод, которым надо решить задачу
    accum -- результат решения задачи на некотором шаге
    st -- шаг интегрирования

    Конструктор для решения задачи численного интегрирования
    MyIntegretionProblem(method::String, ic::Real, st::Real) 

        method -- метод, которым надо решить задачу
        ic -- начальное условие
        st -- шаг интегрирования

"""
const ODEMethodsForIntegretion = Dict("Forward Euler" => :ForwardEuler, "Backward Euler" => :BackwardEuler, 
    "Trapezoidal" => :Trapezoidal)

mutable struct MyIntegretionProblem{Method} <: AbstractMyProblem
    method::String
    accum::Float64
    st::Float64

    function MyIntegretionProblem(method::String, ic::Real, st::Real) 
        haskey(ODEMethodsForIntegretion, method) || 
            throw(ArgumentError("Неподдерживаемый метод интегрирования уравнений"))
        st <= 0 && throw()

        new{ODEMethodsForIntegretion[method]}(method, ic, st)
    end
end