"""
Solve u_t = u_xx in (0,1) with dirichlet BC (u = 0)
using FTCS scheme
"""

import numpy as np
import matplotlib.pyplot as plt

xmin = 0.0; xmax = 1.0

def uexact(x,t):
    return np.exp(-1.0*np.pi**2*t) * np.sin(np.pi*x)

def solve(Tf, cfl, N):
    dx = (xmax - xmin)/(N-1)
    dt = dx**2 * cfl
    iteration = 0

    x = np.linspace(xmin, xmax, N)
    u = uexact(x,0.0); t = 0.0
    fig = plt.figure()
    ax = fig.add_subplot()
    line, = ax.plot(x,u, label="Num sol")        # plot initial condition and assign plot object to "line" variable
    line1, = ax.plot(x,uexact(x,0.0), label="Exact")        # exact sol plotting
    plt.title("Solution")
    plt.xlabel("x"); plt.ylabel("u")
    plt.grid(True); plt.legend()
    plt.draw(); plt.pause(0.1)


    while t < Tf:
        if t + dt > Tf:
            dt = Tf - t
        u[1:-1] = cfl*u[0:-2] + (1-2.0*cfl)*u[1:-1] + cfl*u[2:]
        t += dt; iteration += 1

        line.set_ydata(u)
        line1.set_ydata(uexact(x,t))       # exact sol plotting for every t
        plt.draw(); plt.pause(0.1)
    return iteration


iteration = solve(0.01, 0.2, 30)
print("Iteration = ", iteration)