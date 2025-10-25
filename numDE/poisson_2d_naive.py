# Instead of using matrix notation, this code has loops for calculating u (It's very slow and not recommended).

# without vectorization
from numpy import sin, pi, linspace, meshgrid, zeros, abs, sqrt, mean, linalg, log
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# rhs
f = lambda x,y: 2.0* pi**2 * sin(pi*x) * sin(pi*y)

# Exact sol
uexact = lambda x,y : sin(pi*x) * sin(pi*y)

# domain
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0

# Mesh size
mesh_sizes = [(4,4), (8,8), (16,16), (32,32), (64,64)]
# nx = 50
# ny = 50
final_mse2 = []

for (nx,ny) in mesh_sizes:
  dx, dy = (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1)
  x = linspace(xmin, xmax, nx)
  y = linspace(ymin, ymax, ny)

  print("Grid size = ", nx, "x", ny)
  print("dx, dy    = ", dx, dy)

  # Initialize solution and set boundary conditions
  u = zeros((nx, ny))

  # uexact calculation
  uexact_ = zeros((nx,ny))
  for j in range(0,ny):
      for i in range(0,nx):
          uexact_[i,j] = uexact(x[i],y[j])

  # Iterative solver
  max_iter = 10000
  tol = 1e-6
  error = tol + 1

  for iter in range(max_iter):
      u_old = u.copy()
      squared_error_sum = 0
      for i in range(1, nx - 1):
          for j in range(1, ny - 1):
              u[i, j] = 0.25 * (u_old[i+1, j] + u_old[i-1, j] + u_old[i, j+1] + u_old[i, j-1] + (dx*dy) * f(x[i], y[j]))
              squared_error_sum += (u[i,j] - u_old[i,j])**2
      error = sqrt(squared_error_sum / (nx-2)*(ny-2))
      if error < tol:
          print(f"Converged after {iter+1} iterations.\n")
          break

  # mean square error
  final_mse2.append(linalg.norm(u-uexact_, ord=2))

h = []
for (nx,ny) in mesh_sizes:
  h.append((xmax-xmin)/(nx-1)) # assuming dx = dy (spacing on both axes)

rate2 = log(final_mse2[-1] / final_mse2[-2]) / log(final_mse2[-2] / final_mse2[-3])
print(rate2)

# convergence plot
plt.loglog(h, final_mse2, '*--', label="W/o vectorization")
plt.title("Convergence rate")
plt.legend()
plt.xlabel("h"); plt.ylabel("Rate")


# Color plot error
plt.figure()
plt.title("Error for 64x64")
X, Y = meshgrid(x,y, indexing='ij')
cs = plt.contourf(X,Y, abs(u-uexact_), levels=30)
plt.colorbar(cs)
plt.xlabel("x"); plt.ylabel("y")


plt.show()