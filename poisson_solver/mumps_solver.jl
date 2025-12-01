using MUMPS, MPI, SparseArrays, LinearAlgebra, Plots

# Solving: -Laplace(u) = 1 in (0,1) × (0,1)
#                     u = 0 on x = 0, ∀ y
#                     ∂u/∂η = 0 on x = 1, ∀ y  (first trying for all dirichlet condition)
#                     u = 0 on y = 0, ∀ x
#                     u = 0 on y = 1, ∀ x

using SparseArrays
using SparseMatricesCSR
using Printf
using Plots.PlotMeasures 

# ------------- MPI Initialization -----------------
MPI.Init()
gr()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

# Mesh size
nx = 200
ny = 200

if (rank == 0)
    # rhs
    function f(x,y)
        return sin(2.0*π*x) * sin(2.0*π*y)    
    end

    # Exact solution
    function uexact_(x,y)
        return 1.0/(2.0*(2.0*π)^2) * sin(2.0*π*x) * sin(2.0*π*y)     
    end

    # domain
    xmin, xmax = 0.0, 1.0
    ymin, ymax = 0.0, 1.0


    dx, dy = (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1)
    x = LinRange(xmin, xmax, nx)
    y = LinRange(ymin, ymax, ny)
    # println(length(x))
    # X, Y = meshgrid(x,y, indexing='ij')

    println("Grid size = ", nx, "x", ny)
    println("dx = ", dx)
    println("dy = ", dy)

    # No. of interior points 
    mx, my = nx-2, ny-2     # in case of all Dirichlet condition

    Dxx = sparse((1.0/dx^2) .* Tridiagonal(ones(mx-1), -2.0*ones(mx), ones(mx-1)))
    Dyy = sparse((1.0/dy^2) .* Tridiagonal(ones(my-1), -2.0*ones(my), ones(my-1)))
    Ix, Iy = I(mx), I(my)

    A = -kron(Iy, Dxx) - kron(Dyy, Ix) # fact from numerical analysis (tensor product notation)
    # display(A[1:6,1:6])
    A_rowmajor = SparseMatrixCSR(A)         # to make it same as Praveen sir's code (row major)
    # irn, jcn, a = findnz(A_rowmajor)
    # println(a)
    # println(irn[1:5])
    # println(jcn[1:5])
    # println(a[1:5])
    # X = [xi for yi in y, xi in x]
    # Y = [yi for yi in y, xi in x]

    # rhs vector

    # b = zeros((mx,my))      # loops can be avoided altogether by using X,Y whenever needed
    # for j in range(2,my)
    #     for i in range(2,mx)
    #         b[i,j] = f(x[i],y[j])
    #     end
    # end

    rhs = zeros(mx*my)
    for j in 2:my+1
        for i in 2:mx+1
            rhs[(j-2)*mx + i-1] = f(x[i],y[j])
        end
    end
    # ----------------------- exact solution -----------------------------
    uexact = zeros((nx,ny))
    for i in 1:nx
        for j in 1:nx
            uexact[i,j] = uexact_(x[i],y[j])
        end
    end
    # uexact = uexact_.(X,Y)
    # println(size(uexact))
end

if (rank !== 0)
    A = zeros(nx-2, ny-2)
    rhs = zeros(nx-2)
end

A = MPI.bcast(A, comm, root = 0)
rhs = MPI.bcast(rhs, comm, root = 0)
println("Broadcasting done from rank 0.")

# b = f.(X,Y)
# println(size(b))
# println(b[1:5])


# ----------------------------- MPI PART -----------------------------

# A = sprand(10, 10, 0.2) + I
# A = [1.0 0.0; 0.0 2.0]
# rhs = rand(10)
# rhs = [1, 4]
_sol = solve(A, rhs)
if (rank == 0)
    sol = zeros(nx, ny)
    sol[2:nx-1, 2:ny-1] .= reshape(_sol, nx-2, ny-2) # column major ordering (default in julia)
    # --------------- Plots.jl ---------------------------
    default(size=(600, 500))
    p1 = contour(x, y, transpose(sol), cmap=:viridis, levels=20, xlabel="x", ylabel="y", title="Solution", linewidths=1.2, right_margin = 10mm)
    A = abs.(transpose(sol .- uexact))
    # println("Press Enter to close..."); readline() 
    default(size=(600, 500))
    p2 = contourf(x, y, A, cmap=:viridis, colorbar=true,levels=20, xlabel="x", ylabel="y", title="Absolute Error", right_margin = 10mm)
    # --------------- Plots.jl end -----------------------
    
    # ----------------------- PyPlots.jl -----------------------------
    # figure(figsize=(7,5), constrained_layout=true)
    # A = abs.(transpose(sol .- uexact))
    # p2 = contourf(x, y, A; levels=30, cmap="viridis", antialiased=true)
    # colorbar(p2)                            # Matplotlib will show offset text automatically
    # p1 = contour(x, y, transpose(sol); levels=20, colors="k", linewidths=0.6)
    # xlabel("x"); ylabel("y"); title("Absolute Error")
    # ----------------------- Pyplot.jl -------------------------------



    savefig(p1, "plot1.png")
    savefig(p2, "plot2.png")

    println(maximum(sol .- uexact))
    println("Press Enter to close..."); readline()
end
# println(norm(x - A \ rhs) / norm(x))
# finalize()
MPI.Finalize()