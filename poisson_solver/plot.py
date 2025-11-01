import matplotlib.pyplot as plt
from numpy import meshgrid, linspace, loadtxt, zeros, sin, pi, pad, savetxt

uexact = lambda x,y : 1.0/(2.0*(2.0*pi)**2) * sin(2.0*pi*x) * sin(2.0*pi*y)

# domain
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0

# Mesh size
nx = 10
ny = 10
mx = nx-2; my = ny-2

dx, dy = (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1)
x = linspace(xmin, xmax, nx)
y = linspace(ymin, ymax, ny)
X, Y = meshgrid(x,y, indexing='ij')

sol = loadtxt("sol.txt")
sol = pad(sol, ((1,1), (1,1)), mode='constant', constant_values=0.0)
savetxt("sol_50x50.txt", sol, fmt="%0.6f")
print(sol.shape)

uexact_ = zeros((nx,ny))    # loops can be avoided altogether by using X,Y whenever needed
for j in range(1,my+1):
    for i in range(1,mx+1):
        uexact_[i,j] = uexact(x[i],y[j])
print('Max error = ', abs(sol-uexact_).max())


# Contour plot solution
plt.figure(figsize=(5,5))
plt.title("Solution")
plt.contour(X,Y,sol,20)
plt.xlabel("x"); plt.ylabel("y")

# Color plot error
plt.figure()
plt.title("Error")
cs = plt.contourf(X,Y, abs(sol-uexact_), levels=30)
plt.colorbar(cs)
plt.xlabel("x"); plt.ylabel("y")

plt.show()
