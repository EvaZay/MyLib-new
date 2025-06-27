
@testset "Solve" begin
    @testset "Integration" begin
        @testset "Case1" begin
            method="Forward Euler"
            init=1
            x=[0, 1, 2, 3, 4]
            f(x)=x+1
            result=0
            pr = MyIntegretionProblem(method, init, 1)
            for i in eachindex(x)
                result=calculate(pr, x[i], f(x[i]))
                recalculate!(pr, x[i], f(x[i]))
            end
            @test isapprox(result, 11)
        end

        @testset "Case2" begin
            method="Backward Euler"
            init=1
            x=[0, 1, 2, 3, 4]
            f(x)=x+1
            result=0
            pr = MyIntegretionProblem(method, init, 1)
            for i in eachindex(x)
                result=calculate(pr, x[i], f(x[i]))
                recalculate!(pr, x[i], f(x[i]))
            end
            @test isapprox(result, 16)
        end

        @testset "Case3" begin
            f(x) = x
            x=0
            init=0
            method="Backward Euler"
            st=100//9000
            n = 9000
            result=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 5000.5555555555556)

            f(x) = x
            x=0
            init=0
            method="Forward Euler"
            st=100//9000
            n = 9000
            result=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 4999.444444444444)
        end

        @testset "Case4" begin
            f(x) = x^2
            method="Forward Euler"
            x = -5
            init=0
            st=20//9000
            n = 9000
            result=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 1166.444460905350; atol=1e-9)
        end

        @testset "Case5" begin
            f(x) = x^5+x+10
            method="Forward Euler"
            x = -10.0
            st = 40//9000
            n = 9000
            result=0.0
            init=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 121279917.60664415)
        end

        @testset "Case6" begin
            f(x) = x^5+x+10
            method="Forward Euler"
            x = -10.0
            st = 10//9000
            n = 9000
            result=0.0
            init=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, -166672.2329217837)
        end

        @testset "Case7" begin
            f(x) = x^5-3*x^3+x^2-50
            method="Trapezoidal"
            x = -5
            st = 15//9000
            n = 9000
            result=0.0
            init=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 156653.94820138885)
        end

        @testset "Case8" begin
            f(x) = x+10
            method="Trapezoidal"
            x = -5
            st = 5//10000
            n = 10000
            result=0.0
            init=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 37.501250)
        end

        @testset "Case9" begin
            f(x) = sin(x)
            method="Trapezoidal"
            x = -5
            st = 15//5000
            n = 5000
            result=0.0
            init=0
            pr = MyIntegretionProblem(method, init, st)
            for i=1:n+1
                result=calculate(pr, x, f(x))
                recalculate!(pr, x, f(x))
                x+=st
            end
            @test isapprox(result, 1.1241712589012622)
        end
    end

    @testset "ODE" begin
        @testset "Case1" begin
            f(x,y) = x^2-2*y
            g(x)=3/(4*exp(2*x))+(x^2)/2-x/2+1/4
            method="Euler"
            x=0
            init=1
            n=100
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-1)
            end
        end

        @testset "Case2" begin
            u(x,y, z) = z-2*x
            v(x,y,z)=z
            g(x)=-3*exp(x-3/2)+x^2+2*x-49/(20)
            funcs=[u, v]
            method="Euler"
            x=1.5
            init=[2, -0.2]
            n=10
            st=0.00000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case3" begin
            f(x,y) = x+y
            g(x)=2*exp(x)-x-1
            method="Euler"
            x=0
            init=1
            n=10
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-1)
            end
        end

        @testset "Case13" begin
            f(x,y) = x+y
            g(x)=2*exp(x)-x-1
            method="Runge-Kutta 2"
            x=0
            init=1
            n=100
            st=0.0000001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case14" begin
            u(x,y, z) = z-2*x
            v(x,y,z)=z
            g(x)=-3*exp(x-3/2)+x^2+2*x-49/(20)
            funcs=[u, v]
            method="Runge-Kutta 2"
            x=1.5
            init=[2, -0.2]
            n=10
            st=0.00000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case15" begin
            u(x,y, z, w, k) = (9*x-3*z)/10
            v(x,y,z, w, k)=z
            q(x,y,z, w, k)=0
            f(x,y,z, w, k)=0
            g(x)=-11000/(27*exp((3*x)/(10)))+(x^4)/8-(5*x^3)/3+(58*x^2)/3-(1064*x)/9+11027/(27)
            funcs=[u, v, q, f]
            method="Runge-Kutta 2"
            x=0
            init=[1, 2, 4, 1]
            n=100
            st=0.00000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case7" begin
            f(x,y) = y
            g(x)=exp(x)
            method="Runge-Kutta 3"
            x=0
            init=1
            n=10
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case9" begin
            f(x,y) = cos(x)
            g(x)=sin(x)+1
            method="Runge-Kutta 3"
            x=0
            init=1
            n=100
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case10" begin
            u(x,y, z) = (cos(x)+z)/2
            v(x,y,z)=z
            g(x)=(sin(x))/5-(2*cos(x))/5+(102*exp(x/2))/5-19
            funcs=[u, v]
            method="Runge-Kutta 3"
            x=0
            init=[10, 1]
            n=100
            st=0.0000000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case4" begin
            f(x,y) = x+y
            g(x)=2*exp(x)-x-1
            method="Runge-Kutta 4"
            x=0
            init=1
            n=10
            st=0.00001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case5" begin
            f(x,y) = x^4+y-10
            g(x)=15*exp(x)-x^4-4*x^3-12*x^2-24*x-14
            method="Runge-Kutta 4"
            x=0
            init=1
            n=10
            st=0.000001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case6" begin
            f(x,y) = y
            g(x)=exp(x)
            method="Runge-Kutta 4"
            x=0
            init=1
            n=100
            st=0.000001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case8" begin
            u(x,y,z,w) = 15*x+z-2*w
            v(x,y,z,w)=z
            k(x,y,z,w)=w
            g(x)=((61*sqrt(2)+88)*exp(sqrt(2)*x-x))/2-((61*sqrt(2)-88)*exp(-sqrt(2)*x-x))/2-(15*x^2)/2-30*x-87
            funcs=[u, v, k]
            method="Runge-Kutta 4"
            x=0
            init=[5,4,1]
            n=1000
            st=0.00000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case11" begin
            u(x,y, z, w, k) = (9*x-3*z)/10
            v(x,y,z, w, k)=z
            q(x,y,z, w, k)=0
            f(x,y,z, w, k)=0
            g(x)=-11000/(27*exp((3*x)/(10)))+(x^4)/8-(5*x^3)/3+(58*x^2)/3-(1064*x)/9+11027/(27)
            funcs=[u, v, q, f]
            method="Runge-Kutta 5"
            x=0
            init=[1, 2, 4, 1]
            n=100
            st=0.00000001
            result=0.0

            pr = MyODEProblem(funcs, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

        @testset "Case12" begin
            f(x,y) = x+y
            g(x)=2*exp(x)-x-1
            method="Runge-Kutta 5"
            x=0
            init=1
            n=10
            st=0.000001
            result=0.0

            pr = MyODEProblem(f, method, x, init, st)
            for i=1:n+1
                analitic_solution=g(x)
                result=calculate(pr, x)
                recalculate!(pr, x)
                x+=st
                @test isapprox(result, analitic_solution, atol=1e-4)
            end
        end

    #     @testset "Case16" begin
    #         f(x,y) = cos(x)
    #         g(x)=sin(x)+1
    #         method="Runge-Kutta 4"
    #         x=0
    #         y=1
    #         n=100
    #         h=0.01
    

    #         pr = MyODEProblem(f, method, x, y, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)
    #         calculator_a_tol!(result, g)
    #         @test all(result.a_tol .<= 1e-4)
    #         calculator_a_tol!(result, analitic_solution)
    #         @test all(result.a_tol .<= 1e-4)
    #     end

    #     @testset "Case17" begin
    #         u(x,y, z) = (cos(x)+z)/2
    #         v(x,y,z)=z
    #         g(x)=(sin(x))/5-(2*cos(x))/5+(102*exp(x/2))/5-19
    #         funcs=[u, v]
    #         method="Runge-Kutta 3"
    #         x=0
    #         init=[10, 1]
    #         n=100
    #         h=0.0000000001
    
    #         pr = MyODEProblem(funcs, method, x, init, n, h)
    #         numeric_solution=solve(pr)
    #         analitic_solution=g.(numeric_solution.points)
    #         calculator_a_tol!(numeric_solution, g)
    #         @test all(numeric_solution.a_tol .<= 1e-4)
    #         calculator_a_tol!(numeric_solution, analitic_solution)
    #         @test all(numeric_solution.a_tol .<= 1e-4)
    #     end

    #     @testset "Case18" begin
    #         f(x,y) = y-x
    #         g(x)=x+1
    #         method="Euler"
    #         x=0
    #         y=1
    #         n=100
    #         h=0.01
    

    #         pr = MyODEProblem(f, method, x, y, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)
    #         calculator_a_tol!(result, g)
    #         @test all(result.a_tol .<= 1e-1)
    #         calculator_a_tol!(result, analitic_solution)
    #         @test all(result.a_tol .<= 1e-1)
    #     end

    #     @testset "Case19" begin
    #         u(x,y, z) = 2*z-4*x+5*y
    #         v(x, y, z)=z
    #         funcs=[u,v]
    #         method="Runge-Kutta 4"
    #         x=0
    #         init=[5,2.0]
    #         n=1000000
    #         h=0.000001

    #         g(x)=((47*sqrt(6)+348)*exp(sqrt(6)*x+x))/(300)-((47*sqrt(6)-348)*exp(x-sqrt(6)*x))/(300)+(4*x)/5-8/(25)
    
    #         pr = MyODEProblem(funcs, method, x, init, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)   
    #         analitic_solution=g.(result.points)
    #         calculator_a_tol!(result, g)
    #         @test all(result.a_tol .<= 1e-4)
    #         calculator_a_tol!(result, analitic_solution)
    #         @test all(result.a_tol .<= 1e-4)

    #         @test_nowarn comparison_of_methods(pr, analitic_solution=g, plot_a_tol=true, savefig_a_tol=true)
    #     end

    #     @testset "Case20" begin
    #         u(x,y, z) = z-2*x
    #         v(x,y,z)=z
    #         funcs=[u, v]
    #         method="Runge-Kutta 4"
    #         x=1.5
    #         init=[2, -0.2]
            
    #         n=10
    #         h=0.00000001

    #         g(x)=-3*exp(x-3/2)+x^2+2*x-49/(20)
    
    #         pr = MyODEProblem(funcs, method, x, init, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)   
    #         for i in eachindex(result.points)
    #             @test isapprox(result.values[i], analitic_solution[i], atol=1e-4)
    #         end
    #     end

    #     @testset "Case21" begin
    #         u(x,y, z) = z-2*x-y
    #         v(x, y, z)=z
    #         funcs=[u,v]
    #         method="Runge-Kutta 4"
    #         x=0
    #         init=[2,-0.2]
    #         n=10
    #         h=0.00000001

    #         g(x)=(31*exp(x/2)*sin((sqrt(3)*x)/2))/(5*sqrt(3))+(9*exp(x/2)*cos((sqrt(3)*x)/2))/5-2*x-2
    
    #         pr = MyODEProblem(funcs, method, x, init, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)   
    #         for i in eachindex(result.points)
    #             @test isapprox(result.values[i], analitic_solution[i], atol=1e-4)
    #         end
    #     end

    #     @testset "Case22" begin
    #         u(x,y, z) = 2*z-4*x+5*y
    #         v(x, y, z)=z
    #         funcs=[u,v]
    #         method="Runge-Kutta 4"
    #         x=0
    #         init=[5,2.0]
    #         n=10
    #         h=0.00000001

    #         g(x)=y=((47*sqrt(6)+348)*exp(sqrt(6)*x+x))/(300)-((47*sqrt(6)-348)*exp(x-sqrt(6)*x))/(300)+(4*x)/5-8/(25)
    
    #         pr = MyODEProblem(funcs, method, x, init, n, h)
    #         result=solve(pr)
    #         analitic_solution=g.(result.points)   
    #         for i in eachindex(result.points)
    #             @test isapprox(result.values[i], analitic_solution[i], atol=1e-4)
    #         end
    #     end

        # @testset "" begin
        #     f(x, y) = x^2+3
        #     g(x) = x^3/3+3x+1

        #     method="Euler"
        #     x=0
        #     y=1
        #     n=100
        #     h=0.01
        #     f1(u,p,t)=-2*t


        #     pr = MyODEProblem(f, method, x, y, n, h)
        #     x1, y1 = solve(pr)

        #     method="Runge-Kutta 4"
        #     pr = MyODEProblem(f, method, x, y, n, h)
        #     _, y2 = solve(pr)

        #     y3 = [g(x1[i]) for i in eachindex(x1)]

        #     pl = plot(x1, [y1 y2], title="Solving differential equations", label=["Euler" "Runge-Kutta 4"],  linewidth=3)
        #     plot!(pl, x1, y3, label="Analitics", lw=3, ls=:dot)
        #     xlabel!("x")
        #     ylabel!("y")
        #     savefig(pl, "Plot.png")

        #     abs_tol1=abs.(y3 .- y1)
        #     abs_tol2=abs.(y3 .- y2)
        #     max1=max(abs_tol1...)
        #     max2=max(abs_tol2...)
        #     println(abs_tol2)
        #     if max1>max2 
        #         println("Точнее Рунге-Кутта 4")
        #     else
        #         println("Точнее Эйлер")
        #     end

        #     rel_tol1=(abs_tol1 ./ y3)*100
        #     rel_tol2=(abs_tol2 ./ y3)*100

        #     pogr1=plot(x1, [abs_tol1 abs_tol2], title="Solving differential equations", label=["Euler" "Runge-Kutta 4"], 
        #         linewidth=3)
        #     ylims!(-1e-8, 1e-3)

        #     pogr2=plot(x1, [rel_tol1 rel_tol2], title="Solving differential equations", label=["Euler" "Runge-Kutta 4"], 
        #         linewidth=3)

        #     savefig(pogr1, "Plot1.png")
        #     savefig(pogr2, "Plot2.png")
        # end
     end
end