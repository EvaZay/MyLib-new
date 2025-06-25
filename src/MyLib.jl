module MyLib
    using Plots
    export calculate, recalculate!, comparison_of_methods, plot_result, calculator_a_tol!, calculator_r_tol!, MyODEProblem, MyODEProblem2, MyIntegretionProblem, MyInterpolationProblem, MyODEResult
    include("my_problems.jl")
    include("my_results.jl")
    include("my_solve.jl")
end
