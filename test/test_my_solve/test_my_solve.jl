
@testset "Solve" begin
    @testset "Integration" begin
        @testset "Case1" begin
            method = "Forward Euler"
            init = 1
            x = [0, 1, 2, 3, 4]
            f(x) = x+1
            result = 0
            pr = MyIntegretionProblem(method, init, 1)
            for i in eachindex(x)
                result=calculate(pr, x[i], f(x[i]))
                recalculate!(pr, x[i], f(x[i]))
            end
            @test isapprox(result, 11)
        end

        @testset "Case2" begin
            method = "Backward Euler"
            init = 1
            x = [0, 1, 2, 3, 4]
            f(x) = x + 1
            result = 0
            pr = MyIntegretionProblem(method, init, 1)
            for i in eachindex(x)
                result=calculate(pr, x[i], f(x[i]))
                recalculate!(pr, x[i], f(x[i]))
            end
            @test isapprox(result, 16)
        end

        @testset "Case3" begin
            f(x) = x
            x = 0
            init = 0
            method = "Backward Euler"
            st = 100//9000
            n = 9000
            result = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result = calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 5000.5555555555556)

            f(x) = x
            x = 0
            init = 0
            method = "Forward Euler"
            st = 100//9000
            n = 9000
            result = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 4999.444444444444)
        end

        @testset "Case4" begin
            f(x) = x^2
            method="Forward Euler"
            x = -5
            init = 0
            st = 20//9000
            n = 9000
            result = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result = calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 1166.444460905350; atol=1e-9)
        end

        @testset "Case5" begin
            f(x) = x^5+x+10
            method="Forward Euler"
            x = -10.0
            st = 40//9000
            n = 9000
            result = 0.0
            init = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result = calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 121279917.60664415)
        end

        @testset "Case6" begin
            f(x) = x^5+x+10
            method="Forward Euler"
            x = -10.0
            st = 10//9000
            n = 9000
            result = 0.0
            init = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result = calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, -166672.2329217837)
        end

        @testset "Case7" begin
            f(x) = x^5-3*x^3+x^2-50
            method="Trapezoidal"
            x = -5
            st = 15//9000
            n = 9000
            result = 0.0
            init = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result = calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 156653.94820138885)
        end

        @testset "Case8" begin
            f(x) = x+10
            method="Trapezoidal"
            x = -5
            st = 5//10000
            n = 10000
            result = 0.0
            init = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 37.501250)
        end

        @testset "Case9" begin
            f(x) = sin(x)
            method="Trapezoidal"
            x = -5
            st = 15//5000
            n = 5000
            result = 0.0
            init = 0
            pr = MyIntegretionProblem(method, init, st)
            for i = 1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x += st
            end
            @test isapprox(result, 1.1241712589012622)
        end
    end

    @testset "ODE" begin
        @testset "Case1" begin
            f(x,y) = x^2-2*y
            g(x) = 3/(4*exp(2*x))+(x^2)/2-x/2+1/4
            method = "Euler"
            x = 0
            init = 1
            n = 1000
            st = 0.01
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-1)
            end
        end

        @testset "Case2" begin
            f(x,y) = -y
            g(x) = 1/(exp(x))
            method = "Euler"
            x = 0
            init = 1
            n = 10000
            st = 0.001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-2)
            end
        end

        @testset "Case2" begin
            u(x,y,z) = z-2*x
            v(x,y,z) = z
            g(x) = x^2+2*x+2
            funcs = [u, v]
            method = "Runge-Kutta 4"
            x = -1
            init = [0,1]
            n = 100000
            st = 0.00001
            result = 0.0

            pr = MyODEProblem(funcs, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case3" begin
            f(x,y) = y-x
            g(x) = x+1
            method = "Euler"
            x = 0
            init = 1
            n = 1000
            st = 0.001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_a_tol!(pr, points, results, g)
            @test all(pr.a_tol .<= 1e-3)
            calculator_a_tol!(pr, results, analitic_solution)
            @test all(pr.a_tol .<= 1e-3)
        end

        @testset "Case4" begin
            f(x,y) = x+y
            g(x) = 2*exp(x)-x-1
            method = "Runge-Kutta 2"
            x = 0
            init = 1
            n = 10000
            st = 0.0001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-2)
            end
        end

        @testset "Case5" begin
            f(x,y) = -y
            g(x) = 1/(exp(x))
            method = "Runge-Kutta 2"
            x = 0
            init = 1
            n = 10000
            st = 0.001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

         @testset "Case6" begin
            f(x,y) = cos(x)
            g(x) = sin(x)+1
            method="Runge-Kutta 2"
            x = 0
            init = 1
            n = 10000
            st = 0.0001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_a_tol!(pr, points, results, g)
            @test all(pr.a_tol .<= 1e-3)
            calculator_a_tol!(pr, results, analitic_solution)
            @test all(pr.a_tol .<= 1e-3)
        end

        @testset "Case7" begin
            f(x,y) = y-x
            g(x) = x+1
            method="Runge-Kutta 2"
            x = 0
            init = 1
            n = 1000
            st = 0.001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_a_tol!(pr, points, results, g)
            @test all(pr.a_tol .<= 1e-3)
            calculator_a_tol!(pr, results, analitic_solution)
            @test all(pr.a_tol .<= 1e-3)
        end

        @testset "Case8" begin
            f(x,y) = y
            g(x)=exp(x)
            method="Runge-Kutta 3"
            x=0
            init=1
            n=100000
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case9" begin
            f(x,y) = -y
            g(x) = 1/(exp(x))
            method = "Runge-Kutta 3"
            x = 0
            init = 1
            n = 10000
            st = 0.001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case10" begin
            f(x,y) = x+y
            g(x) = 2*exp(x)-x-1
            method="Runge-Kutta 3"
            x = 0
            init = 1
            n = 10000
            st = 0.0001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case11" begin
            f(x,y) = y-x
            g(x) = x+1
            method = "Runge-Kutta 3"
            x = 0
            init = 1
            n = 1000
            st = 0.001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_a_tol!(pr, points, results, g)
            @test all(pr.a_tol .<= 1e-3)
            calculator_a_tol!(pr, results, analitic_solution)
            @test all(pr.a_tol .<= 1e-3)
        end

        @testset "Case12" begin
            f(x,y) = -y
            g(x) = 1/(exp(x))
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 10000
            st = 0.001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case13" begin
            f(x,y) = x+y
            g(x) = 2*exp(x)-x-1
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 10000
            st = 0.0001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x +=st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case14" begin
            f(x,y) = x^4+y-10
            g(x) = 15*exp(x)-x^4-4*x^3-12*x^2-24*x-14
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 100000
            st = 0.00001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case15" begin
            f(x,y) = y
            g(x) = exp(x)
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 10000
            st = 0.0001
            result = 0.0

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                analitic_solution = g(x)
                result = calculate(pr, x)
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution, atol=1e-3)
            end
        end

        @testset "Case16" begin
            f(x,y) = y-x
            g(x) = x+1
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 1000
            st = 0.001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_a_tol!(pr, points, results, g)
            @test all(pr.a_tol .<= 1e-3)
            calculator_a_tol!(pr, results, analitic_solution)
            @test all(pr.a_tol .<= 1e-3)
        end

        @testset "Case17" begin
            f(x,y) = y-x
            g(x) = x+1
            method = "Runge-Kutta 4"
            x = 0
            init = 1
            n = 1000
            st = 0.001
            result = 0.0
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(f, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_r_tol!(pr, points, results, g)
            @test all(pr.r_tol .<= 1e-3)
            calculator_r_tol!(pr, results, analitic_solution)
            @test all(pr.r_tol .<= 1e-3)
        end

        @testset "Case18" begin
            u(x,y,z) = z-2*x
            v(x,y,z) = z
            g(x) = x^2+2*x+2
            funcs = [u, v]
            x = -1
            init = [0,1]
            n = 100000
            st = 0.00001

            result = 0.0
            method = "Runge-Kutta 3"
            analitic_solution = zeros(Float64, n + 1)
            points = zeros(Float64, n + 1)
            results = zeros(Float64, n + 1)

            pr = MyODEProblem(funcs, method, x, init, st, n)
            for i = 1:n+1
                points[i] = x
                analitic_solution[i] = g(x)
                result = calculate(pr, x)
                results[i] = result
                recalculate!(pr, x)
                x += st
                @test isapprox(result, analitic_solution[i], atol=1e-3)
            end
            calculator_r_tol!(pr, points, results, g)
            @test all(pr.r_tol .<= 1e-3)
            calculator_r_tol!(pr, results, analitic_solution)
            @test all(pr.r_tol .<= 1e-3)
        end

        @testset "Case19" begin
            u(x,y, z) = 2*z-4*x+5*y
            v(x, y, z)=z
            funcs=[u,v]
            method="Runge-Kutta 4"
            x=0
            init=[2, 5]
            n=1000000
            st=0.000001

            pr = MyODEProblem(funcs, method, x, init, st, n)

            g(x) = -((89*sqrt(6)-774)*exp(sqrt(6)*x+x))/(300)+((89*sqrt(6)+774)*exp(x-sqrt(6)*x))/(300)+(2*x)/5-4/(25)
            @test_nowarn comparison_of_methods(funcs, x, init, st, n, plot_results = false, savefig_results = false)

            results = zeros(Float64, n)
            points = zeros(Float64, n)
            for i = 1:n
                results[i] = g(x)
                points[i] = x
                x += st
            end
            @test_nowarn plot_result(pr, points, results)
        end

        @testset "Case20" begin
            f(x,y) = -y
            g(x) = 1/(exp(x))
            x = 0
            init = 1
            n = 10000
            st = 0.001
            result = 0.0
            method = "Euler"

            pr = MyODEProblem(f, method, x, init, st, n)

            @test_nowarn comparison_of_methods(f, x, init, st, n, plot_results = false, savefig_results = false)

            results = zeros(Float64, n)
            points = zeros(Float64, n)
            for i = 1:n
                results[i] = g(x)
                points[i] = x
                x += st
            end
            @test_nowarn plot_result(pr, points, results)
        end

        @testset "Case21" begin
            u(x,y,z) = z-2*x
            v(x,y,z) = z
            g(x) = x^2+2*x+2
            funcs = [u, v]
            x = -1
            init = [0,1]
            n = 100000
            st = 0.00001
            result = 0.0
            method = "Runge-Kutta 4"

            pr = MyODEProblem(funcs, method, x, init, st, n)

            @test_nowarn comparison_of_methods(funcs, x, init, st, n, plot_results = false, savefig_results = false)

            results = zeros(Float64, n)
            points = zeros(Float64, n)
            for i = 1:n
                results[i] = g(x)
                points[i] = x
                x += st
            end
            @test_nowarn plot_result(pr, points, results)
        end
     end
end