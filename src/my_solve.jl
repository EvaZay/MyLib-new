function recalculate!(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    pr.accum+=pr.st*y
end

function calculate(pr::MyIntegretionProblem{:BackwardEuler}, x, y)
    return pr.accum+pr.st*y
end


function recalculate!(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    pr.accum += pr.st * y
end

function calculate(pr::MyIntegretionProblem{:ForwardEuler}, x, y)
    return pr.accum
end


function recalculate!(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    pr.accum += pr.st * y
end

function calculate(pr::MyIntegretionProblem{:Trapezoidal}, x, y)
    return pr.accum + pr.st*y / 2
end


function recalculate!(pr::MyODEProblem{:Euler}, x)
    N=length(pr.results)
    for i=1:N
        pr.results[i]+=+pr.st*(pr.coeffs1[i])
    end
end

function calculate(pr::MyODEProblem{:Euler}, x)
    N=length(pr.results)
    for i=1:N
        pr.coeffs1[i]=pr.funcs[i](x, pr.results...)
    end

    for j=1:N
        result=pr.results[j]+pr.st*(pr.coeffs1[j])

        if j==N
            return result
        end
    end
end


function recalculate!(pr::MyODEProblem{:RungeKutta2}, x)
    N=length(pr.results)
    for i=1:N
        pr.results[i]+=pr.st/2*(pr.coeffs1[i]+pr.coeffs2[i])
    end
end

function calculate(pr::MyODEProblem{:RungeKutta2}, x)
    N=length(pr.results)
    for i=1:N
        pr.coeffs1[i]=pr.funcs[i](x, pr.results...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs1[end-k+1] for k in 1:N])
        pr.coeffs2[i]=pr.funcs[i](x+pr.st, args...)
    end

    for j=1:N
        result=pr.results[j]+pr.st/2*(pr.coeffs1[j]+pr.coeffs2[j])

        if j==N
            return result
        end
    end
end


function recalculate!(pr::MyODEProblem{:RungeKutta3}, x)
    N=length(pr.results)
    for i=1:N
        pr.results[i]+=pr.st/6*(pr.coeffs1[i]+4*pr.coeffs2[i]+pr.coeffs3[i])
    end
end

function calculate(pr::MyODEProblem{:RungeKutta3}, x)
    N=length(pr.results)
    for i=1:N
        pr.coeffs1[i]=pr.funcs[i](x, pr.results...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs1[end-k+1]*pr.st/2 for k in 1:N])
        pr.coeffs2[i]=pr.funcs[i](x+pr.st/2, args...)
    end

    for i=1:N
        args = Tuple([pr.results[k]-pr.coeffs1[end-k+1]+2*pr.coeffs2[end-k+1]*pr.st/2 for k in 1:N])
        pr.coeffs3[i]=pr.funcs[i](x+pr.st, args...)
    end

    for j=1:N
        result=pr.results[j]+pr.st/6*(pr.coeffs1[j]+4*pr.coeffs2[j]+pr.coeffs3[j])

        if j==N
            return result
        end
    end
end


function recalculate!(pr::MyODEProblem{:RungeKutta4}, x)
    N=length(pr.results)
    for i=1:N
        pr.results[i]+=pr.st/6*(pr.coeffs1[i]+2*pr.coeffs2[i]+2*pr.coeffs3[i]+pr.coeffs4[i])
    end
end

function calculate(pr::MyODEProblem{:RungeKutta4}, x)
    N=length(pr.results)
    for i=1:N
        pr.coeffs1[i]=pr.funcs[i](x, pr.results...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs1[end-k+1]*pr.st/2 for k in 1:N])
        pr.coeffs2[i]=pr.funcs[i](x+pr.st/2, args...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs2[end-k+1]*pr.st/2 for k in 1:N])
        pr.coeffs3[i]=pr.funcs[i](x+pr.st/2, args...)
    end

    for i=1:N
        args = Tuple([pr.results[k]-pr.coeffs1[end-k+1]+pr.coeffs3[end-k+1]*pr.st for k in 1:N])
        pr.coeffs4[i]=pr.funcs[i](x+pr.st, args...)
    end

    for j=1:N
        result=pr.results[j]+pr.st/6*(pr.coeffs1[j]+2*pr.coeffs2[j]+2*pr.coeffs3[j]+pr.coeffs4[j])

        if j==N
            return result
        end
    end
end


function recalculate!(pr::MyODEProblem{:RungeKutta5}, x)
    N=length(pr.results)
    for i=1:N
        pr.results[i]+=pr.st/6*(pr.coeffs1[i]+4*pr.coeffs4[i]+pr.coeffs5[i])
    end
end

function calculate(pr::MyODEProblem{:RungeKutta5}, x)
    N=length(pr.results)
    for i=1:N
        pr.coeffs1[i]=pr.funcs[i](x, pr.results...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs1[end-k+1]*pr.st for k in 1:N])
        pr.coeffs2[i]=pr.funcs[i](x+pr.st/3, args...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs2[end-k+1]*pr.st/2+pr.coeffs1[end-k+1]*pr.st/2 for k in 1:N])
        pr.coeffs3[i]=pr.funcs[i](x+pr.st/3, args...)
    end

    for i=1:N
        args = Tuple([pr.results[k]+pr.coeffs3[end-k+1]*9/8*pr.st+2/8*pr.coeffs1[end-k+1]*pr.st for k in 1:N])
        pr.coeffs4[i]=pr.funcs[i](x+pr.st/2, args...)
    end
    for i=1:N
        args = Tuple([pr.results[k]-9/2*pr.coeffs3[end-k+1]*pr.st+6*pr.coeffs4[end-k+1]*pr.st+3/2*pr.coeffs1[end-k+1]*pr.st for k in 1:N])
        pr.coeffs4[i]=pr.funcs[i](x+pr.st, args...)
    end

    for j=1:N
        result=pr.results[j]+pr.st/6*(pr.coeffs1[j]+4*pr.coeffs4[j]+pr.coeffs5[j])

        if j==N
            return result
        end
    end
end

# function calculator_a_tol!(
#     result::MyODEResult,
#     analitic_solution::Function
#     )
#     analitic_values = [analitic_solution(result.points[i]) for i in eachindex(result.points)]
#     result.a_tol.=abs.(analitic_values .- result.values)
# end
# function calculator_a_tol!(
#     result::MyODEResult,
#     analitic_solution::Vector{<:Real}
#     )
#     result.a_tol.=abs.(analitic_solution .- result.values)
# end

# function calculator_r_tol!(
#     result::MyODEResult,
#     analitic_solution::Function
#     )
#     if all(isnan.result.a_tol) 
#         calculator_a_tol(result, analitic_solution)
#     end
#     analitic_values = [analitic_solution(result.points[i]) for i in eachindex(result.points)]
#     result.r_tol.=result.a_tol ./ analitic_values
# end
# function calculator_r_tol!(
#     result::MyODEResult,
#     analitic_solution::Vector{<:Real}
#     )
#     if all(isnan.result.a_tol) 
#         calculator_a_tol(result, analitic_solution)
#     end
#     result.r_tol.=result.a_tol ./ analitic_values
# end

# function comparison_of_methods(
#     pr::MyODEProblem;
#     methods::Vector{String}=["Euler","Runge-Kutta 2","Runge-Kutta 3","Runge-Kutta 4","Runge-Kutta 5"],
#     plot_results::Bool=true,
#     savefig_results::Bool=true,
#     filename_result::String="plot",
#     analitic_solution::Union{Function, Nothing}=nothing,
#     plot_a_tol::Bool=false,
#     savefig_a_tol::Bool=false,
#     filename_a_tol::String="a_tol",
#     plot_r_tol::Bool=false,
#     savefig_r_tol::Bool=false,
#     filename_r_tol::String="r_tol"
# )
#     (plot_a_tol || plot_r_tol) && isnothing(analitic_solution) && throw(ArgumentError("Когда параметр plot_a_tol или plot_r_tol имеет значение true, \
#     должно быть передано аналитическое решение"))
#     results=MyODEResult[]
#     pl_res = plot() 
#     xlabel!("x")
#     ylabel!("y")
#     pl_a_tol=plot()
#     xlabel!("x")
#     ylabel!("a_tol")
#     pl_r_tol=plot()
#     xlabel!("x")
#     ylabel!("r_tol")
 
#     for method in methods

#         prob = MyODEProblem(pr.funcs, method, pr.x[1], [result[1] for result in pr.results], pr.n, pr.h)
#         res=solve(prob)

#         if plot_a_tol
#             calculator_a_tol!(res, analitic_solution)
#         end

#         if plot_r_tol
#             calculator_r_tol!(res, analitic_solution)
#         end

#         if plot_results 
#             plot!(pl_res, res.points, res.values, label=method, lw=3, ls=:dot)
#             if savefig_results
#                 savefig(pl_res, filename_result*".png")
#             end
#         end

#         if plot_a_tol
#             plot!(pl_a_tol, res.points, res.a_tol, label=method, lw=3, ls=:dot)
#             if savefig_a_tol
#                 savefig(pl_a_tol, filename_a_tol*".png")
#             end
#         end

#         if plot_r_tol
#             plot!(pl_r_tol, res.points, res.r_tol, label=method, lw=3, ls=:dot)
#             if savefig_r_tol
#                 savefig(pl_r_tol, filename_r_tol*".png")
#             end
#         end
#     end
# end

# function plot_result(
#     result::MyODEResult; 
#     plot_result::Bool=true, 
#     plot_a_tol::Bool=false, 
#     plot_r_tol::Bool=false,
#     savefig_result::Bool=false,
#     savefig_a_tol::Bool=false,
#     savefig_r_tol::Bool=false,
#     filename_result::String="plot",
#     filename_a_tol::String="a_tol",
#     filename_r_tol::String="r_tol"
#     )

#     if plot_result
#         pl_res = plot(result.points, result.values, title="Solving differential equation", linewidth=3)
#         xlabel!("x")
#         ylabel!("y")

#         if savefig_result
#             savefig(pl_res, filename_result*".png")
#         end
#     end

#     if plot_a_tol
#         pl_a_tol = plot(result.points, result.a_tol, title="Absolute tolerance", linewidth=3)
#         xlabel!("x")
#         ylabel!("a_tol")

#         if savefig_a_tol
#             savefig(pl_a_tol, filename_a_tol*".png")
#         end
#     end

#     if plot_r_tol
#         pl_r_tol = plot(result.points, result.r_tol, title="Relative tolerance", linewidth=3)
#         xlabel!("x")
#         ylabel!("r_tol")

#         if savefig_r_tol
#             savefig(pl_r_tol, filename_r_tol*".png")
#         end
#     end
# end