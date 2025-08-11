"""
Подсчет результата для задачи численного интегрирования 
calculate(pr::MyIntegretionProblem, x, y)

pr -- Экземпляр задачи
x -- координата x
y -- значение функции
"""

function calculate(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    return pr.accum + pr.st*y
end

function calculate(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    return pr.accum
end

function calculate(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    return pr.accum + pr.st*y / 2
end


"""
Пересчет внутренних параметров для задачи численного интегрирования 
recalculate!(pr::MyIntegretionProblem, x, y)

pr -- Экземпляр задачи
x -- координата x
y -- значение функции
"""
function recalculate!(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    pr.accum += pr.st * y
end

function recalculate!(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    pr.accum += pr.st * y
end

function recalculate!(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    pr.accum += pr.st * y
end


"""
Подсчет результата для задачи численного решения дифференциальных уравнений n-ого порядка
calculate(pr::MyODEProblem, x)

pr -- Экземпляр задачи
x -- координата x
"""
function calculate(pr::MyODEProblem{:Euler}, x)
    N = length(pr.results)
    for i = 1:N
        args = Tuple([pr.results[end-k+1] for k in 1:N])
        pr.coeffs1[i] = pr.st * pr.funcs[i](x, args...)
    end

    for j = 1:N
        result = pr.results[j] + pr.coeffs1[j]

        if j == N
            return result
        end
    end
end

function calculate(pr::MyODEProblem{:RungeKutta2}, x)
    N = length(pr.results)
    for i = 1:N
        args = Tuple([pr.results[end-k+1] for k in 1:N])
        pr.coeffs1[i] = pr.st * pr.funcs[i](x, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] + pr.coeffs1[k] for k in 1:N])
        pr.coeffs2[i] = pr.st * pr.funcs[i](x + pr.st, args...)
    end

    for j = 1:N
        result = pr.results[j] + (pr.coeffs1[j] + pr.coeffs2[j]) / 2

        if j == N
            return result
        end
    end
end

function calculate(pr::MyODEProblem{:RungeKutta3}, x)
    N = length(pr.results)
    for i = 1:N
        args = Tuple([pr.results[end-k+1] for k in 1:N])
        pr.coeffs1[i] = pr.st * pr.funcs[i](x, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] + pr.coeffs1[k] / 2 for k in 1:N])
        pr.coeffs2[i] = pr.st * pr.funcs[i](x + pr.st / 2, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] - pr.coeffs1[k] + 2 * pr.coeffs2[k] for k in 1:N])
        pr.coeffs3[i] = pr.st * pr.funcs[i](x + pr.st, args...)
    end

    for j = 1:N
        result = pr.results[j] +  (pr.coeffs1[j] + 4 * pr.coeffs2[j] + pr.coeffs3[j]) / 6

        if j == N
            return result
        end
    end
end

function calculate(pr::MyODEProblem{:RungeKutta4}, x)
    N = length(pr.results)
    for i = 1:N
        args = Tuple([pr.results[end-k+1] for k in 1:N])
        pr.coeffs1[i] = pr.st * pr.funcs[i](x, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] + pr.coeffs1[k] / 2 for k in 1:N])
        pr.coeffs2[i]=pr.st * pr.funcs[i](x + pr.st / 2, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] + pr.coeffs2[k] / 2 for k in 1:N])
        pr.coeffs3[i] = pr.st * pr.funcs[i](x + pr.st / 2, args...)
    end

    for i = 1:N
        args = Tuple([pr.results[end-k+1] + pr.coeffs3[k] for k in 1:N])
        pr.coeffs4[i] = pr.st * pr.funcs[i](x + pr.st, args...)
    end

    for j = 1:N
        result = pr.results[j] + (pr.coeffs1[j] + 2 * pr.coeffs2[j] + 2 * pr.coeffs3[j] + pr.coeffs4[j])/6

        if j == N
            return result
        end
    end
end

"""
Пересчет внутренних параметров для задачи численного решения дифференциальных уравнений n-ого порядка
recalculate!(pr::MyODEProblem, x)

pr -- Экземпляр задачи
x -- координата x
"""
function recalculate!(pr::MyODEProblem{:Euler}, x)
    N = length(pr.results)
    for i = 1:N
        pr.results[i] += (pr.coeffs1[i])
    end
end

function recalculate!(pr::MyODEProblem{:RungeKutta2}, x)
    N = length(pr.results)
    for i = 1:N
        pr.results[i] += (pr.coeffs1[i] + pr.coeffs2[i]) / 2
    end
end

function recalculate!(pr::MyODEProblem{:RungeKutta3}, x)
    N = length(pr.results)
    for i = 1:N
        pr.results[i] += (pr.coeffs1[i] + 4 * pr.coeffs2[i] + pr.coeffs3[i]) / 6
    end
end

function recalculate!(pr::MyODEProblem{:RungeKutta4}, x)
    N = length(pr.results)
    for i = 1:N
        pr.results[i] += (pr.coeffs1[i] + 2 * pr.coeffs2[i] + 2 * pr.coeffs3[i] + pr.coeffs4[i]) / 6
    end
end

"""
Метод подсчета абсолютной погрешности
function calculator_a_tol!(pr, x, result, analitic_solution)

    pr -- задача
    x -- значение x на некотором шаге
    result -- решение задачи на некотром шаге
    analitic_solution -- аналитическое решение задачи в виде функции

function calculator_a_tol!(pr, result, analitic_solution)

    pr -- задача
    result -- решение задачи на некотром шаге
    analitic_solution -- аналитическое решение задачи в виде массива значений

"""

function calculator_a_tol!(
    pr::MyODEProblem,
    x::Vector{<:Real},
    result::Vector{<:Real}, 
    analitic_solution::Function
    )
    analitic_values = [analitic_solution(x[i]) for i in eachindex(x)]
    pr.a_tol .= abs.(analitic_values .- result)
end
function calculator_a_tol!(
    pr::MyODEProblem,
    result::Vector{<:Real},
    analitic_solution::Vector{<:Real}
    )
    pr.a_tol .= abs.(analitic_solution .- result)
end

"""
Метод подсчета относительной погрешности
function calculator_r_tol!(pr, x, result, analitic_solution)

    pr -- задача
    x -- значение x на некотором шаге
    result -- решение задачи на некотром шаге
    analitic_solution -- аналитическое решение задачи в виде функции

function calculator_r_tol!(pr, result, analitic_solution)

    pr -- задача
    result -- решение задачи на некотром шаге
    analitic_solution -- аналитическое решение задачи в виде массива значений

"""

function calculator_r_tol!(
    pr::MyODEProblem,
    x::Vector{<:Real},
    result::Vector{<:Real}, 
    analitic_solution::Function
    )
    #Добавить тест, который будет заходить в этот иф
    if all(isnan, pr.a_tol) 
        calculator_a_tol!(pr, x, result, analitic_solution)
    end
    analitic_values = [analitic_solution(x[i]) for i in eachindex(x)]
    pr.r_tol .= pr.a_tol ./ analitic_values
end
function calculator_r_tol!(
    pr::MyODEProblem,
    result::Vector{<:Real},
    analitic_solution::Vector{<:Real}
    )
    #Добавить тест, который будет заходить в этот иф
    if all(isnan, pr.a_tol) 
        calculator_a_tol!(pr, result, analitic_solution)
    end
    pr.r_tol .= pr.a_tol ./ analitic_solution
end

"""
Метод нахождения решений задачи введенными методами, построение графиков решений, абсолютной и относительной погрешностей

function comparison_of_methods(
    funcs::Vector{<:Function}, x, init::Vector{<:Real}, st, n; methods, plot_results, savefig_results, filename_result, analitic_solution, plot_a_tol, savefig_a_tol, 
    filename_a_tol, plot_r_ol, savefig_r_tol, filename_r_tol
)
    funcs -- система ОДУ
    x -- начальное значение x
    init -- начальные условия
    st -- шаг дифференцирования
    n -- количество шагов
    methods -- вектор методов, которыми необходимо решить задачу
    plot_results -- необходимость в построении графика решений задачи (по умолчанию стоит true)
    savefig_results -- необходимость в сохранении графика решений задачи (по умолчанию стоит true)
    filename_result -- имя файла с графиком решений задачи (по умолчанию стоит "plot")
    analitic_solution -- аналитическое решение задачи в виде функции
    plot_a_tol -- необходимость в построении графика абсолютных погрешностей задачи (по умолчанию стоит false)
    savefig_a_tol -- необходимость в сохранении графика абсолютных погрешностей задачи (по умолчанию стоит false)
    filename_a_tol -- имя файла с графиком абсолютных погрешностей задачи (по умолчанию стоит "a_tol")
    plot_r_tol -- необходимость в построении графика относительных погрешностей задачи (по умолчанию стоит false)
    savefig_r_tol -- необходимость в сохранении графика относительных погрешностей задачи (по умолчанию стоит false)
    filename_r_tol -- имя файла с графиком относительных погрешностей задачи (по умолчанию стоит "r_tol")
"""

function comparison_of_methods(
    funcs::Vector{<:Function},
    x::Real,
    init::Vector{<:Real},
    st::Real,
    n::Int64;
    methods::Vector{String}=["Euler","Runge-Kutta 2","Runge-Kutta 3","Runge-Kutta 4"],
    plot_results::Bool=true,
    savefig_results::Bool=true,
    filename_result::String="plot",
    plot_analitic_solution::Bool=false,
    analitic_solution::Union{Function, Nothing}=nothing,
    plot_a_tol::Bool=false,
    savefig_a_tol::Bool=false,
    filename_a_tol::String="a_tol",
    plot_r_tol::Bool=false,
    savefig_r_tol::Bool=false,
    filename_r_tol::String="r_tol"
)
    (plot_a_tol || plot_r_tol) && isnothing(analitic_solution) && 
        throw(ArgumentError("Когда параметр plot_a_tol или plot_r_tol имеет значение true, \
            должно быть передано аналитическое решение"))
    pl_res = plot() 
    xlabel!("x")
    ylabel!("y")
    pl_a_tol = plot()
    xlabel!("x")
    ylabel!("a_tol")
    pl_r_tol = plot()
    xlabel!("x")
    ylabel!("r_tol")

    analitic_sol = zeros(Float64, n + 1)
    points = zeros(Float64, n + 1)
 
    for method in methods
        results = zeros(Float64, n + 1)
        point = x

        prob = MyODEProblem(funcs, method, x, init, st, n)
        for i = 1:n+1
            points[i] = point
            if analitic_solution != nothing
                analitic_sol[i] = analitic_solution(point)
            end
            results[i] = calculate(prob, point)
            recalculate!(prob, point)
            point += st
        end

        if plot_a_tol
            calculator_a_tol!(prob, results, analitic_sol)
        end

        if plot_r_tol
            calculator_r_tol!(prob, results, analitic_sol)
        end

        if plot_results 
            plot!(pl_res, points, results, label = method, lw = 3, ls = :dot)
            if savefig_results
                savefig(pl_res, filename_result * ".png")
            end
        end

        if plot_a_tol
            plot!(pl_a_tol, points, prob.a_tol, label = method, lw = 3, ls = :dot)
            if savefig_a_tol
                savefig(pl_a_tol, filename_a_tol * ".png")
            end
        end

        if plot_r_tol
            plot!(pl_r_tol, points, prob.r_tol, label = method, lw = 3, ls = :dot)
            if savefig_r_tol
                savefig(pl_r_tol, filename_r_tol * ".png")
            end
        end
    end

    if plot_analitic_solution
        point = x
        for i = 1:n+1
            points[i] = point
            analitic_sol[i] = analitic_solution(point)
            point += st
        end
        plot!(pl_res, points, analitic_sol, label = "analitic solution", lw = 3, ls = :dot)
        savefig(pl_res, filename_result * ".png")
    end
    
end

"""
Метод нахождения решений задачи введенными методами, построение графиков решений, 
абсолютной и относительной погрешностей

function comparison_of_methods(
    funcs::Function, x, init::Real, st, n; methods, plot_results, 
    savefig_results, filename_result, analitic_solution, plot_a_tol, savefig_a_tol, 
    filename_a_tol, plot_r_ol, savefig_r_tol, filename_r_tol
)
    funcs -- функция, которую необходимо продефференцировать
    x -- начальное значение x
    init -- начальное условие
    st -- шаг дифференцирования
    n -- количество шагов
    methods -- вектор методов, которыми необходимо решить задачу
    plot_results -- необходимость в построении графика решений задачи (по умолчанию стоит true)
    savefig_results -- необходимость в сохранении графика решений задачи (по умолчанию стоит true)
    filename_result -- имя файла с графиком решений задачи (по умолчанию стоит "plot")
    analitic_solution -- аналитическое решение задачи в виде функции
    plot_a_tol -- необходимость в построении графика абсолютных погрешностей задачи (по умолчанию стоит false)
    savefig_a_tol -- необходимость в сохранении графика абсолютных погрешностей задачи (по умолчанию стоит false)
    filename_a_tol -- имя файла с графиком абсолютных погрешностей задачи (по умолчанию стоит "a_tol")
    plot_r_tol -- необходимость в построении графика относительных погрешностей задачи (по умолчанию стоит false)
    savefig_r_tol -- необходимость в сохранении графика относительных погрешностей задачи (по умолчанию стоит false)
    filename_r_tol -- имя файла с графиком относительных погрешностей задачи (по умолчанию стоит "r_tol")
"""

function comparison_of_methods(
    funcs::Function,
    x::Real,
    init::Real,
    st::Real,
    n::Int64;
    methods::Vector{String}=["Euler","Runge-Kutta 2","Runge-Kutta 3","Runge-Kutta 4"],
    plot_results::Bool=true,
    savefig_results::Bool=true,
    filename_result::String="plot",
    plot_analitic_solution::Bool=false,
    analitic_solution::Union{Function, Nothing}=nothing,
    plot_a_tol::Bool=false,
    savefig_a_tol::Bool=false,
    filename_a_tol::String="a_tol",
    plot_r_tol::Bool=false,
    savefig_r_tol::Bool=false,
    filename_r_tol::String="r_tol"
)
    (plot_a_tol || plot_r_tol) && isnothing(analitic_solution) && 
        throw(ArgumentError("Когда параметр plot_a_tol или plot_r_tol имеет значение true, \
            должно быть передано аналитическое решение"))
    pl_res = plot() 
    xlabel!("x")
    ylabel!("y")
    pl_a_tol = plot()
    xlabel!("x")
    ylabel!("a_tol")
    pl_r_tol = plot()
    xlabel!("x")
    ylabel!("r_tol")
    analitic_sol = zeros(Float64, n + 1)
    points = zeros(Float64, n + 1)
    for method in methods
        results=zeros(Float64, n + 1)
        point = x

        prob = MyODEProblem(funcs, method, x, init, st, n)
        for i = 1:n+1
            points[i] = point
            if analitic_solution != nothing
                analitic_sol[i] = analitic_solution(point)
            end
            results[i] = calculate(prob, point)
            recalculate!(prob, point)
            point += st
        end

        if plot_a_tol
            calculator_a_tol!(prob, results, analitic_sol)
        end

        if plot_r_tol
            calculator_r_tol!(prob, results, analitic_sol)
        end

        if plot_results 
            plot!(pl_res, points, results, label = method, lw = 3, ls = :dot)
            if savefig_results
                savefig(pl_res, filename_result * ".png")
            end
        end

        if plot_a_tol
            plot!(pl_a_tol, points, prob.a_tol, label = method, lw = 3, ls = :dot)
            if savefig_a_tol
                savefig(pl_a_tol, filename_a_tol * ".png")
            end
        end

        if plot_r_tol
            plot!(pl_r_tol, points, prob.r_tol, label = method, lw = 3, ls = :dot)
            if savefig_r_tol
                savefig(pl_r_tol, filename_r_tol * ".png")
            end
        end
    end

    if plot_analitic_solution 
        point = x
        for i = 1:n+1
            points[i] = point
            analitic_sol[i] = analitic_solution(point)
            point += st
        end
        plot!(pl_res, points, analitic_sol, label = "analitic solution", lw = 3, ls = :dot)
        savefig(pl_res, filename_result * ".png")
    end
end

function plot_result(
    pr::MyODEProblem,
    points::Vector{<:Real},
    results::Vector{<:Real}; 
    plot_result::Bool=true, 
    plot_a_tol::Bool=false, 
    plot_r_tol::Bool=false,
    savefig_result::Bool=true,
    savefig_a_tol::Bool=false,
    savefig_r_tol::Bool=false,
    filename_result::String="plot",
    filename_a_tol::String="a_tol",
    filename_r_tol::String="r_tol"
    )

    if plot_result
        pl_res = plot(points, results, title = "Solving differential equation", lw = 3, ls = :dot)
        xlabel!("x")
        ylabel!("y")

        if savefig_result
            savefig(pl_res, filename_result*".png")
        end
    end

    if plot_a_tol
        pl_a_tol = plot(points, pr.a_tol, title = "Absolute tolerance", linewidth = 3)
        xlabel!("x")
        ylabel!("a_tol")

        if savefig_a_tol
            savefig(pl_a_tol, filename_a_tol * ".png")
        end
    end

    if plot_r_tol
        pl_r_tol = plot(points, pr.r_tol, title = "Relative tolerance", linewidth = 3)
        xlabel!("x")
        ylabel!("r_tol")

        if savefig_r_tol
            savefig(pl_r_tol, filename_r_tol * ".png")
        end
    end
end