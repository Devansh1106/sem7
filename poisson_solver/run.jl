using LinearAlgebra

include("poisson_mumps.jl")
A, b = generate_poisson(5)

# # Flatten A and b into a single stream
# n = size(A, 1)
# data = [Float64(n); vec(A); b]
data = [irn; ]

# Convert to raw binary bytes
bytes = reinterpret(UInt8, data)

# Step 2. Run the MPI solver and pipe data
cmd = `mpiexec -n 4 ./mumps_solver`
io = IOBuffer()

# Launch process with stdin/stdout piping
open(cmd, "r+", stdout=io) do proc
    write(proc, bytes)    # send A and b to solver stdin
    close(proc.in)        # important: signal EOF to C side
end

# Step 3. Read solver output (as Float64 array)
seekstart(io)
result_bytes = read(io)
x = reinterpret(Float64, result_bytes)

println("Solution vector from MPI solver:")
println(x)
