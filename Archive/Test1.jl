using Plots

A = 2;
B = 3;

C = A + B;

A = 1:1:10;

B = [1 2; 3 4];
C = sin(B)

# display(C)

C = sin.(B)

# display(C)


display(plot(A,A))



B = range(0,5,100);


# display("Symulacja gotowaa")
# # p2 = plot(title = "Tytuł")
# for a = range(1,5,1000)
#     D = sin.(B*a)
#     p2= plot(B,D)
#     display(p2)
#     sleep(0)
# end

# sleep(6)


# function  Test(a)
#     a = a+1;
# end

# Test(1);

# using BenchmarkTools

# x = rand(10_000_000)  # Duży wektor

# @btime sum($x)
# @btime dot(ones(length($x)), $x)
# @btime ones(length($x))' * $x


# typeof(H1)
# typeof(H2)
# Iterations.Iteration(NonlinearityFull, Sets.Set1, Sets.VIni1)

# A = 3;

# Plot(LaplaceDiscretisation.SCB')




# trace = scatter(x=x, y=y, mode="lines")
# layout = Layout(title="Wykres sin(x*j)")
# fig = plot(trace, layout)
# display(fig)


# global running = true

# function test()
#     x = 0.0:0.1:10
#     y = sin.(x)
#     fig = plot(x,y)
#     display(fig)

#     j = 1;
#     # while running
#     while (true)
#         y = sin.(x.*j)
#         f = plot(x,y)
#         display(f)
#         j = j+0.1
#         # linia[1] = y

#         # sleep(0.1)
#         # println(j)
#         sleep(0.01)
        
#     end
# end



# function listen_for_commands()
#     while running
#         print("\n> ")
#         input = readline(stdin)
#         sleep(1)
#         if input == "exit"
#             println("Kończenie...")
#             global running = false
#         else
#             println("Nieznana komenda.")
#         end
#     end
#     print("Koniec")
# end



# Startujemy nasłuchiwanie w osobnym wątku
# Threads.@spawn listen_for_commands
# test()


A = rand(100,100)
B = rand(100)

function T1(n,A,B)
    for i =1:n
        B = A*B;
    end
end

function M(A::Matrix{Float64},B::Vector{Float64})
    return A*B;
end

function T2(n::Int64,A::Matrix{Float64},B::Vector{Float64})
    for i =1:n
        B = A*B;
    end
end

@time T1(100000,A,B)

@time T2(100000,A,B)

using Revise
using LinearAlgebra
using SparseArrays
using Plots
using Statistics
using  BenchmarkTools

includet("ParametersHydra1D.jl")

# Solvers File
includet("SolverSteps.jl")

using  .LaplaceDiscretisation

Mat = Lap.Per8;

N = 10000; # Number of discretization points
dt = 0.0001; # Time step
L = 1; # Domain size
dx = L/N; # Spatial step
x = range(0,L,N); # Discretisation Vector
SicCosNodes = 3; # Sine Cosine Nodes
T = 100.0;

κ = 10.0;
D = 0.1;

V = (200 .*rand(N));


for t = 0.0:dt:T
    t1 = time()
    for i = 1:1000
    V = V + dt.*(D .*Mat*V ./ (dx.^2) + - V + κ * exp.(V) ./ mean(exp.(V)));  
    end
    display(time()-t1)
    display(plot(V))
    sleep(0.01);
end

using LinearSolve
using ..LaplaceDiscretisation

# A = rand(100,100);
A = 100 .*Lap.Neu8 - I(100);
B = rand(100);

function Brute()
    return A^(-1)*B;

end
@time Brute();

C = A^(-1);

function Brute2(C)
    return C*B;

end

@time Brute2(C);

function Medium()
    solve(LinearProblem(A,B));
end

@time Medium();


AS = sparse(A);

function Sparse()
    solve(LinearProblem(AS,B));
end

@time Sparse();

function select_matrix(Matrix, Lap)
    cases = Dict(
        "Dirichlet 2"  => Lap.Dir2,
        "Dirichlet 4"  => Lap.Dir4,
        "Dirichlet 6"  => Lap.Dir6,
        "Dirichlet 8"  => Lap.Dir8,
        "Neumann 2"    => Lap.Neu2,
        "Neumann 4"    => Lap.Neu4,
        "Neumann 6"    => Lap.Neu6,
        "Neumann 8"    => Lap.Neu8,
        "Periodic 2"   => Lap.Per2,
        "Periodic 4"   => Lap.Per4,
        "Periodic 6"   => Lap.Per6,
        "Periodic 8"   => Lap.Per8
    )
    
    return get(cases, Matrix, error("Nieznana wartość Matrix: $Matrix"))
end


using PyCall, BenchmarkTools

# Rozmiar macierzy
N = 2000;

# Tworzenie losowych macierzy w Julii
A = rand(N, N);
B = rand(N, N);

function MullJulia(A::Matrix{Float64}, B::Matrix{Float64})
        A = A * B;
    return A;
end

function MullPython(A::PyObject, B::PyObject)
        A = np.matmul(A, B);
    return A;
end

A1 = PyObject(A);
B1 = PyObject(B);

# Pomiar czasu mnożenia macierzy w Julii
println("Czas wykonania w Julii:")
@btime MullJulia(A,B);

println("Czas wykonania w Pythonie (NumPy):")
@btime MullPython(A1,B1);


using FEniCS
mesh = UnitSquareMesh(8,8)
V = FunctionSpace(mesh,"P",1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]", degree=2)
u = TrialFunction(V)
bc1 = DirichletBC(V,u_D, "on_boundary")
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u),grad(v))*dx
L = f*v*dx
U = FeFunction(V)
lvsolve(a,L,U,bc1) #linear variational solver
errornorm(u_D, U, norm="L2")
get_array(L) #this returns an array for the stiffness matrix
get_array(U) #this returns an array for the solution values
vtkfile = File("poisson/solution.pvd")
vtkfile << U.pyobject #exports the solution to a vtkfile