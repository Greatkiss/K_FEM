import numpy as np
#the physical property
global rho
global mu
global nu
rho = 998
mu = 1
nu = mu/rho
dt = 0.1
dx = 0.1

def estimated_u(Nx, Ny, u, v, p) :
    ucap = u
    for i in range(1, Nx) :
        for j in range(1, Ny) :
            advection = -dt/dx*(u[i,j]*(u[i+1,j] - u[i,j]) + v[i,j]*(u[i, j+1] - u[i,j]))
            viscosity = nu*dt*(1/dx**2*(u[i+1,j] - 2*u[i,j] + u[i-1,j]) + 1/dx**2*(u[i,j+1] - 2*u[i,j] + u[i,j-1]))
            pressure =  -dt/dx*rho*(p[i+1,j] - p[i,j])
            ucap[i,j] = u[i,j] + advection + viscosity + pressure
    return ucap

def estimated_v(Nx, Ny, u, v, p) :
    vcap = v
    for i in range(1, Nx) :
        for j in range(1, Ny) :
            advection = -dt/dx*(u[i,j]*(v[i+1,j] - v[i,j]) + v[i,j]*(v[i, j+1] - v[i,j]))
            viscosity = nu*dt*(1/dx**2*(v[i+1,j] - 2*v[i,j] + v[i-1,j]) + 1/dx**2*(v[i,j+1] - 2*v[i,j] + v[i,j-1]))
            pressure =  -dt/dx*rho*(p[i,j+1] - p[i,j])
            vcap[i,j] = v[i,j] + advection + viscosity + pressure
    return vcap

def sol_pdash(Nx, Ny, u, v):
    pdash = np.zeros((Nx+2, Ny+2))
    pddash = np.zeros((Nx+2, Ny+2))
    for i in range(1, Nx) :
        for j in range(1, Ny) :
            f = rho/dt*( (u[i+1,j] - u[i,j])/dx + (v[i,j+1] - v[i,j])/dx )
            pddash[i,j] = (pdash[i+1,j] + pdash[i-1,j] + pdash[i,j+1] + pdash[i,j-1] -f*dx**2)/4
    zazansa = 1
    while zazansa > 0.001:
        zansa = 0
        pdash = pddash
        print(pdash)
        for i in range(1, Nx) :
            for j in range(1, Ny) :
                f = rho/dt*( (u[i+1,j] - u[i,j])/dx + (v[i,j+1] - v[i,j])/dx )
                pddash[i,j] = (pdash[i+1,j] + pdash[i-1,j] + pdash[i,j+1] + pdash[i,j-1] -f*dx**2)/4
                zansa = zansa + abs((pddash[i,j] - pdash[i,j]))
        zazansa = zansa
        print(pddash)
    return pddash 

def sol_udash(Nx, Ny, pdash):
    udash = np.zeros((Nx+2, Ny+2))
    for i in range(1, Nx) :
        for j in range(1, Ny) :
            udash[i,j] = -dt/rho*(pdash[i+1,j] - pdash[i,j])/dx
    return udash

def sol_vdash(Nx, Ny, pdash):
    vdash = np.zeros((Nx+2, Ny+2))
    for i in range(1, Nx) :
        for j in range(1, Ny) :
            vdash[i,j] = -dt/rho*(pdash[i,j+1] - pdash[i,j])/dx
    return vdash