# Inspired from repo: https://github.com/cpraveen/numpde/blob/master/bvp2d/bvp2da.py 

"""
Poisson solving: -Laplace(u) = f
                    u = 0
"""
from numpy import sin, pi, linspace, meshgrid, zeros, ones, reshape, abs
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# rhs
f = lambda x,y: sin(2.0*pi*x) * sin(2.0*pi*y)

# Exact sol
uexact = lambda x,y : 1.0/(2.0*(2.0*pi)**2) * sin(2.0*pi*x) * sin(2.0*pi*y)

# domain
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0

# Mesh size
nx = 50
ny = 50

dx, dy = (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1)
x = linspace(xmin, xmax, nx)
y = linspace(ymin, ymax, ny)
X, Y = meshgrid(x,y, indexing='ij')

print("Grid size = ", nx, "x", ny)
print("dx, dy    = ", dx, dy)

# No. of interior points
mx, my = nx-2, ny-2

Dxx = (1.0/dx**2) * diags([ones(mx), -2.0*ones(mx), ones(mx)],[-1,0,1],(mx,mx))
Dyy = (1.0/dy**2) * diags([ones(my), -2.0*ones(my), ones(my)], [-1,0,1],(my,my))
Ix, Iy = eye(mx), eye(my)

A = -kron(Iy, Dxx) - kron(Dyy, Ix) # fact from numerical analysis (tensor product notation)

# rhs vector
b = zeros((mx,my))      # loops can be avoided altogether by using X,Y whenever needed
for j in range(1,my+1):
    for i in range(1,mx+1):
        b[i-1,j-1] = f(x[i],y[j])

b = reshape(b, mx*my, order='F')

# uexact calculation
uexact_ = zeros((nx,ny))    # loops can be avoided altogether by using X,Y whenever needed
for j in range(1,my+1):
    for i in range(1,mx+1):
        uexact_[i,j] = f(x[i],y[j])

# uexact_ = reshape(uexact_, mx*my, order='F')

# solve
sol = spsolve(A,b)

# Reshape to array
u = zeros((nx,ny))  # BCs are 0
u[1:-1,1:-1] = reshape(sol, (mx,my), order='F')
print('Max error = ', abs(u-uexact_).max())


# Contour plot solution
plt.figure(figsize=(5,5))
plt.title("Solution")
plt.contour(X,Y,u,20)
plt.xlabel("x"); plt.ylabel("y")

# Color plot error
plt.figure()
plt.title("Error")
cs = plt.contourf(X,Y, abs(u-uexact_), levels=30)
plt.colorbar(cs)
plt.xlabel("x"); plt.ylabel("y")

plt.show()