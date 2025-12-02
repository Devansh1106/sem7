# Created by: Devansh Tripathi (github.com/Devansh1106)
# Created at: UTC+05:30 23:23 

# Solving: -Laplace(u) = sin(2.0*π*x) * sin(2.0*π*y)  in (0,1) × (0,1)
#                     u = 0 on x = 0, ∀ y
#                     u = 0 on x = 1, ∀ y
#                     u = 0 on y = 0, ∀ x
#                     u = 0 on y = 1, ∀ x
# Exact sol: 1.0/(2.0*(2.0*π)^2) * sin(2.0*π*x) * sin(2.0*π*y)
using MUMPS, MPI, SparseArrays, LinearAlgebra, PyPlot

# ------------- MPI Initialization -----------------
MPI.Init()
# gr()
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
    A_rowmajor = sparse(A)         # to make it same as Praveen sir's code (row major)

    # rhs vector
    rhs = zeros(mx*my)
    @inbounds for j in 2:my+1
        for i in 2:mx+1
            rhs[(j-2)*mx + i-1] = f(x[i],y[j])
        end
    end
    # ----------------------- exact solution -----------------------------
    uexact = zeros((nx,ny))
    @inbounds for i in 1:nx
        for j in 1:nx
            uexact[i,j] = uexact_(x[i],y[j])
        end
    end
end

if (rank !== 0)
    A = zeros(nx-2, ny-2)
    rhs = zeros(nx-2)
end

@time begin 
    A = MPI.bcast(A, comm, root = 0)
    rhs = MPI.bcast(rhs, comm, root = 0)
    println("Broadcasting done from rank 0.")
    _sol = solve(A, rhs)
end

if (rank == 0)
    sol = zeros(nx, ny)
    sol[2:nx-1, 2:ny-1] .= reshape(_sol, nx-2, ny-2) # column major ordering (default in julia)
    # --------------- Plots.jl ---------------------------
    # default(size=(600, 500))
    # p1 = contour(x, y, sol', cmap=:viridis, levels=20, xlabel="x", ylabel="y", title="Solution", linewidths=1.2, right_margin = 10mm)
    # A = abs.(transpose(sol .- uexact))
    # # println(A[50,50])
    # # println("Press Enter to close..."); readline() 
    # p2 = contourf(x, y, A, cmap=:viridis, colorbar=true,levels=30, xlabel="x", ylabel="y", title="Absolute Error", right_margin = 10mm)
    # --------------- Plots.jl end -----------------------

    # ----------------------- PyPlots.jl -----------------------------
    figure(figsize=(6,5), constrained_layout=true)
    p1 = contour(x, y, transpose(sol), cmap="viridis", levels=20, linewidths=0.6)
    savefig("plot1_jl.png", dpi=300)
    Z = abs.(transpose(sol - uexact))
    p2 = contourf(x, y, Z, levels=30, cmap="viridis")
    colorbar(p2)                            # Matplotlib will show offset text automatically
    xlabel("x"); ylabel("y"); title("Absolute Error")
    savefig("plot2_jl.png", dpi=300)
    # ----------------------- Pyplot.jl -------------------------------

    # savefig(p1, "plot1_jl.png")
    # savefig(p2, "plot2_jl.png")
    print("Absolute Error: ")
    println(maximum(sol .- uexact))
    # println("Press Enter to close..."); readline()
end
MPI.Finalize()