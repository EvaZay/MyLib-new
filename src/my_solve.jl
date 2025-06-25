function recalculate!(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    pr.accum+=pr.st*y
end

function calculate(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    return pr.accum+pr.st*y
end



# function solve(pr::MyIntegretionProblem{:BackwardEuler})
#     sum=0.0
#     h=(pr.b-pr.a)/pr.n
#     for i=1:pr.n
#         x=pr.a+(i-1)*h
#         sum+=pr.f(x)
#     end
#     result=h*sum
# end
function recalculate!(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    pr.accum += pr.st * y
end

function calculate(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    return pr.accum
end

# function solve(pr::MyIntegretionProblem{:ForwardEuler})
#     sum=0.0
#     h=(pr.b-pr.a)/pr.n
#     for i=1:pr.n
#         x=pr.a+i*h
#         sum+=pr.f(x)
#     end
#     result=h*sum
# end

function recalculate!(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    pr.accum += pr.st * y
end

function calculate(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    return pr.accum + pr.st*y / 2
end

function solve(pr::MyIntegretionProblem{:Trapezoidal})
    h=(pr.b-pr.a)/pr.n
    sum=0.0
    res1=(pr.f(pr.a)+pr.f(pr.b))/2
    for i=1:pr.n-1
        x=pr.a+i*h
        sum+=pr.f(x)
    end
    result=h*(res1+sum)
end

function solve(pr::MyIntegretionProblem{:Simpsons})
    h=(pr.b-pr.a)/pr.n
    sum=0.0
    for i=1:pr.n
        x_prev=pr.a+(i-1)*h
        x=pr.a+i*h
        sum+=pr.f(x)+4*pr.f((x_prev+x)/2)+pr.f(x_prev)
    end
    result=h*sum/6
end


function solve(pr::MyODEProblem{:Euler})
    n=length(pr.x)-1

    N=length(pr.results)
    
    for i=1:n
        pr.x[i+1]=pr.x[i]+pr.h

        for j=1:N
            args = Tuple([pr.results[k][i] for k in 1:N])
            pr.coeffs1[j]=pr.funcs[j](pr.x[i], args...)
        end
        

        for j=1:N
            pr.results[j][i+1]=pr.results[j][i]+pr.h*(pr.coeffs1[j])
        end
    end

    result=MyODEResult(pr.x, pr.results[end])

    return result
end

function solve(pr::MyODEProblem{:RungeKutta2})
    n=length(pr.x)-1

    N=length(pr.results)
    
    for i=1:n
        pr.x[i+1]=pr.x[i]+pr.h

        for j=1:N
            args = Tuple([pr.results[k][i] for k in 1:N])
            pr.coeffs1[j]=pr.funcs[j](pr.x[i], args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs1[end-k+1] for k in 1:N])
            pr.coeffs2[j]=pr.funcs[j](pr.x[i]+pr.h, args...)
        end


        for j=1:N
            pr.results[j][i+1]=pr.results[j][i]+pr.h/2*(pr.coeffs1[j]+pr.coeffs2[j])
        end
    end

    result=MyODEResult(pr.x, pr.results[end])

    return result
end

function solve(pr::MyODEProblem{:RungeKutta3})
    n=length(pr.x)-1

    N=length(pr.results)
    
    for i=1:n
        pr.x[i+1]=pr.x[i]+pr.h

        for j=1:N
            args = Tuple([pr.results[k][i] for k in 1:N])
            pr.coeffs1[j]=pr.funcs[j](pr.x[i], args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs1[end-k+1]*pr.h/2 for k in 1:N])
            pr.coeffs2[j]=pr.funcs[j](pr.x[i]+pr.h/2, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]-pr.coeffs1[end-k+1]+2*pr.coeffs2[end-k+1]*pr.h/2 for k in 1:N])
            pr.coeffs3[j]=pr.funcs[j](pr.x[i]+pr.h, args...)
        end
        

        for j=1:N
            pr.results[j][i+1]=pr.results[j][i]+pr.h/6*(pr.coeffs1[j]+4*pr.coeffs2[j]+pr.coeffs3[j])
        end
    end

    result=MyODEResult(pr.x, pr.results[end])

    return result
end

function solve(pr::MyODEProblem{:RungeKutta4})
    n=length(pr.x)-1

    N=length(pr.results)
    
    for i=1:n

        pr.x[i+1]=pr.x[i]+pr.h

        for j=1:N
            args = Tuple([pr.results[k][i] for k in 1:N])
            pr.coeffs1[j]=pr.funcs[j](pr.x[i], args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs1[end-k+1]*pr.h/2 for k in 1:N])
            pr.coeffs2[j]=pr.funcs[j](pr.x[i]+pr.h/2, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs2[end-k+1]*pr.h/2 for k in 1:N])
            pr.coeffs3[j]=pr.funcs[j](pr.x[i]+pr.h/2, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs3[end-k+1]*pr.h for k in 1:N])
            pr.coeffs4[j]=pr.funcs[j](pr.x[i]+pr.h, args...)
        end


        for j=1:N
            pr.results[j][i+1]=pr.results[j][i]+pr.h/6*(pr.coeffs1[j]+2*pr.coeffs2[j]+2*pr.coeffs3[j]+pr.coeffs4[j])
        end
    end

    result=MyODEResult(pr.x, pr.results[end])

    return result
end

function solve(pr::MyODEProblem{:RungeKutta5})
    n=length(pr.x)-1

    N=length(pr.results)
    
    for i=1:n

        pr.x[i+1]=pr.x[i]+pr.h

        for j=1:N
            args = Tuple([pr.results[k][i] for k in 1:N])
            pr.coeffs1[j]=pr.funcs[j](pr.x[i], args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs1[end-k+1]*pr.h for k in 1:N])
            pr.coeffs2[j]=pr.funcs[j](pr.x[i]+pr.h/3, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs2[end-k+1]*pr.h/2+pr.coeffs1[end-k+1]*pr.h/2 for k in 1:N])
            pr.coeffs3[j]=pr.funcs[j](pr.x[i]+pr.h/3, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]+pr.coeffs3[end-k+1]*9/8*pr.h+2/8*pr.coeffs1[end-k+1]*pr.h for k in 1:N])
            pr.coeffs4[j]=pr.funcs[j](pr.x[i]+pr.h/2, args...)
        end

        for j=1:N
            args = Tuple([pr.results[k][i]-9/2*pr.coeffs3[end-k+1]*pr.h+6*pr.coeffs4[end-k+1]*pr.h+3/2*pr.coeffs1[end-k+1]*pr.h for k in 1:N])
            pr.coeffs4[j]=pr.funcs[j](pr.x[i]+pr.h, args...)
        end


        for j=1:N
            pr.results[j][i+1]=pr.results[j][i]+pr.h/6*(pr.coeffs1[j]+4*pr.coeffs4[j]+pr.coeffs5[j])
        end
    end

    result=MyODEResult(pr.x, pr.results[end])

    return result
end

function calculator_a_tol!(
    result::MyODEResult,
    analitic_solution::Function
    )
    analitic_values = [analitic_solution(result.points[i]) for i in eachindex(result.points)]
    result.a_tol.=abs.(analitic_values .- result.values)
end
function calculator_a_tol!(
    result::MyODEResult,
    analitic_solution::Vector{<:Real}
    )
    result.a_tol.=abs.(analitic_solution .- result.values)
end

function calculator_r_tol!(
    result::MyODEResult,
    analitic_solution::Function
    )
    if all(isnan.result.a_tol) 
        calculator_a_tol(result, analitic_solution)
    end
    analitic_values = [analitic_solution(result.points[i]) for i in eachindex(result.points)]
    result.r_tol.=result.a_tol ./ analitic_values
end
function calculator_r_tol!(
    result::MyODEResult,
    analitic_solution::Vector{<:Real}
    )
    if all(isnan.result.a_tol) 
        calculator_a_tol(result, analitic_solution)
    end
    result.r_tol.=result.a_tol ./ analitic_values
end

function comparison_of_methods(
    pr::MyODEProblem;
    methods::Vector{String}=["Euler","Runge-Kutta 2","Runge-Kutta 3","Runge-Kutta 4","Runge-Kutta 5"],
    plot_results::Bool=true,
    savefig_results::Bool=true,
    filename_result::String="plot",
    analitic_solution::Union{Function, Nothing}=nothing,
    plot_a_tol::Bool=false,
    savefig_a_tol::Bool=false,
    filename_a_tol::String="a_tol",
    plot_r_tol::Bool=false,
    savefig_r_tol::Bool=false,
    filename_r_tol::String="r_tol"
)
    (plot_a_tol || plot_r_tol) && isnothing(analitic_solution) && throw(ArgumentError("Когда параметр plot_a_tol или plot_r_tol имеет значение true, \
    должно быть передано аналитическое решение"))
    results=MyODEResult[]
    pl_res = plot() 
    xlabel!("x")
    ylabel!("y")
    pl_a_tol=plot()
    xlabel!("x")
    ylabel!("a_tol")
    pl_r_tol=plot()
    xlabel!("x")
    ylabel!("r_tol")
 
    for method in methods

        prob = MyODEProblem(pr.funcs, method, pr.x[1], [result[1] for result in pr.results], pr.n, pr.h)
        res=solve(prob)

        if plot_a_tol
            calculator_a_tol!(res, analitic_solution)
        end

        if plot_r_tol
            calculator_r_tol!(res, analitic_solution)
        end

        if plot_results 
            plot!(pl_res, res.points, res.values, label=method, lw=3, ls=:dot)
            if savefig_results
                savefig(pl_res, filename_result*".png")
            end
        end

        if plot_a_tol
            plot!(pl_a_tol, res.points, res.a_tol, label=method, lw=3, ls=:dot)
            if savefig_a_tol
                savefig(pl_a_tol, filename_a_tol*".png")
            end
        end

        if plot_r_tol
            plot!(pl_r_tol, res.points, res.r_tol, label=method, lw=3, ls=:dot)
            if savefig_r_tol
                savefig(pl_r_tol, filename_r_tol*".png")
            end
        end
    end
end

function plot_result(
    result::MyODEResult; 
    plot_result::Bool=true, 
    plot_a_tol::Bool=false, 
    plot_r_tol::Bool=false,
    savefig_result::Bool=false,
    savefig_a_tol::Bool=false,
    savefig_r_tol::Bool=false,
    filename_result::String="plot",
    filename_a_tol::String="a_tol",
    filename_r_tol::String="r_tol"
    )

    if plot_result
        pl_res = plot(result.points, result.values, title="Solving differential equation", linewidth=3)
        xlabel!("x")
        ylabel!("y")

        if savefig_result
            savefig(pl_res, filename_result*".png")
        end
    end

    if plot_a_tol
        pl_a_tol = plot(result.points, result.a_tol, title="Absolute tolerance", linewidth=3)
        xlabel!("x")
        ylabel!("a_tol")

        if savefig_a_tol
            savefig(pl_a_tol, filename_a_tol*".png")
        end
    end

    if plot_r_tol
        pl_r_tol = plot(result.points, result.r_tol, title="Relative tolerance", linewidth=3)
        xlabel!("x")
        ylabel!("r_tol")

        if savefig_r_tol
            savefig(pl_r_tol, filename_r_tol*".png")
        end
    end
end