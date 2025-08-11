module MyLib

using Plots

export calculate, 
    recalculate!, 
    comparison_of_methods, 
    plot_result, 
    calculator_a_tol!, 
    calculator_r_tol!, 
    MyODEProblem, 
    MyIntegretionProblem


include("my_problems.jl")
include("my_solve.jl")

end
