

# Solving: -Laplace(u) = 1 in (0,1) × (0,1)
#                     u = 0 on x = 0, ∀ y
#                     ∂u/∂η = 0 on x = 1, ∀ y  (first trying for all dirichlet condition)
#                     u = 0 on y = 0, ∀ x
#                     u = 0 on y = 1, ∀ x

using LinearAlgebra
using SparseArrays
using Printf

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

# Mesh size
nx = 6
ny = 6

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
irn, jcn, a = findnz(A)
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

b = zeros(mx*my)
for j in 2:my+1
    for i in 2:mx+1
        b[(j-2)*mx + i-1] = f(x[i],y[j])
    end
end


# b = f.(X,Y)
println(size(b))
println(b[1:5])

# exact solution
uexact = zeros((nx,ny))
for i in 1:nx
    for j in 1:nx
        uexact[i,j] = uexact_(x[i],y[j])
    end
end
# uexact = uexact_.(X,Y)
# println(size(uexact))

nnz = length(a)
# println(nnz)
# Open file and write data
open("matrix.txt", "w") do io
    # First line: n and nnz
    println(io, "$mx $nnz")

    # Matrix entries: i j a
    for k in 1:nnz
        println(io, "$(irn[k]) $(jcn[k]) $(a[k])")
    end
    for k in 1:mx*my
        print(io, "$(b[k]) ")
    end
    # RHS values
    # for i in 1:mx
    #     print(io, b[i])
    #     if i < mx
    #         print(io, " ")
    #     # else
    #     #     println(io)
    #     end
    # end
end

open("exact_sol.txt", "w") do io
#     println(io, "$uexact")
# end
    for i in 1:nx
        for j in 1:nx
            # print(io, "$(uexact[i,j])")
            @printf(io, "%.6f", uexact[i,j])
            print(io, " ")
        end
        println(io)
    end
end
# for k in 1:nnz
#     print(io, "$(rhs[k])")
# end

println("Matrix assembly done!")















# sol_, time_taken = ccall((:say_hi, "./mumps_solver.so"), Ptr{Cdouble}, Cdouble,   # return type
#                                           (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}, # argument types
#                                           Cint, # length of each array)
#                                           Cint), # order of matrix
#                                           irn, jcn, a, b, length(a), mx)
# sol = unsafe_wrap(Vector{Cdouble}, sol_, length(a))

# println(sol[1:5])

# sol, time_taken = run(`mpicc -O2 mumps_solver.c -o mumps_solver`)
# println(time_taken)
# println(sol[1:5])

