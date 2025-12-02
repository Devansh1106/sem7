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
icntl = default_icntl[:]
icntl[1] = -1; icntl[2] = 0; icntl[3] = 0; icntl[4] = 0; # To supress all the MUMPS error messages
mumps = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)    #mumps_symmetric can cause blowing up absolute error values if matrix is actually not symmetric

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

    # rhs
    rhs = zeros(mx*my)
    @inbounds for j in 2:my+1
        for i in 2:mx+1
            rhs[(j-2)*mx + i-1] = f(x[i],y[j])
        end
    end
    # exact solution
    uexact = zeros((nx,ny))
    @inbounds for i in 1:nx
        for j in 1:nx
            uexact[i,j] = uexact_(x[i],y[j])
        end
    end
    associate_matrix!(mumps, A_rowmajor)
    associate_rhs!(mumps, rhs)
end

# main solver part that call MUMPS internally
@time begin
    factorize!(mumps)
    solve!(mumps)
end
MPI.Barrier(comm)

if (rank == 0)
    sol = zeros(nx, ny)
    _sol = get_solution(mumps)
    sol[2:nx-1, 2:ny-1] .= reshape(_sol, nx-2, ny-2) # column major ordering (default in julia)
    # ----------------------- PyPlots.jl -----------------------------
    figure(figsize=(6,5), constrained_layout=true)
    p1 = contour(x, y, transpose(sol), cmap="viridis", levels=20, linewidths=0.6)
    savefig("plot1.png", dpi=300)
    Z = abs.(transpose(sol - uexact))
    p2 = contourf(x, y, Z, levels=30, cmap="viridis")
    colorbar(p2)                            # Matplotlib will show offset text automatically
    savefig("plot2.png", dpi=300)
    xlabel("x"); ylabel("y"); title("Absolute Error")
    # ----------------------- Pyplot.jl end -------------------------------

    println(maximum(abs.(sol - uexact)))
end
finalize(mumps)     # destroy the mumps struct internally
MPI.Finalize()