from ctypes import sizeof
from curses import A_ALTCHARSET
import numpy as np
import math
import matplotlib.pyplot as plt
from mods.mkmesh import mkmesh
from mods.simplemethod import estimated_u, estimated_v, sol_pdash, sol_udash, sol_vdash
from mods.BC import BC

dx = 0.1

#the shape of simulation area
rect = [1,1]
size = [dx, dx]
Nw = int(rect[0] / size[0])
Nh = int(rect[1] / size[1])
u = mkmesh(rect[0], rect[1], size[0], size[1])
v = mkmesh(rect[0], rect[1], size[0], size[1])
p = mkmesh(rect[0], rect[1], size[0], size[1])

# 散布図として描画
fig, ax = plt.subplots(figsize=(10, 10))
x = np.linspace(0, 11, 12)
y = np.linspace(0, 11, 12)
xx, yy = np.meshgrid(x,y)
q = plt.quiver(x, y ,u, v,color='black',angles='xy',scale_units='xy', scale=4.5)
# x, y ラベル
ax.set_xlabel('$X$', fontsize=15)
ax.set_ylabel('$Y$', rotation=0, fontsize=15)

#the physical property
global rho
global mu
global nu
rho = 998
mu = 1
nu = rho/mu

#the B.C.
BC(u, v, p, Nw, Nh)

q = plt.quiver(x, y ,u.T, v.T,color='black',angles='xy',scale_units='xy', scale=4.5)
plt.show()

#main
zansa = 1

for time in range (5):
    #B.C. again
    BC(u, v, p, Nw, Nh)
    udash = np.zeros((Nw+2, Nh+2))
    vdash = np.zeros((Nw+2, Nh+2))
    pdash = np.zeros((Nw+2, Nh+2))
    ucap = estimated_u(Nw, Nh, u, v, p)
    vcap = estimated_v(Nw, Nh, u, v, p)
    pdash = sol_pdash(Nw, Nh, ucap, vcap)
    udash = sol_udash(Nw, Nh, pdash)
    vdash = sol_vdash(Nw, Nh, pdash)
    q = plt.quiver(x, y ,udash.T, vdash.T,color='black',angles='xy',scale_units='xy', scale=4.5)
    plt.show()
    p = p + pdash
    u = u + udash
    v = v + vdash
    zansa = np.sum(pdash) + np.sum(udash)+ np.sum(vdash)
