import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

dt = 0.01
maxt = 1
g = 9.81
nsteps = int(maxt/dt)
r = np.array([0., 0., 0.])
v = np.array([0., 3., 3.])
o = np.array([0., 0., 0.])

rs = np.empty((nsteps, 3))
vs = np.empty((nsteps, 3))
os = np.empty((nsteps, 3))

def Fr(v):
    return v
def Fv(v, o):
    F = np.array([0., 0., -g])
    return F
def Fo(o):
    return np.array([0., 0., 0.])

def SimRK4(r, v, o):
    kr1 = dt*Fr(v)
    kv1 = dt*Fv(v, o)
    ko1 = dt*Fo(o)

    kr2 = dt*Fr(v + kv1/2)
    kv2 = dt*Fv(v + kv1/2, o + ko1/2)
    ko2 = dt*Fo(o + ko1/2)

    kr3 = dt*Fr(v + kv2/2)
    kv3 = dt*Fv(v + kv2/2, o + ko2/2)
    ko3 = dt*Fo(o + ko2/2)

    kr4 = dt*Fr(v + kv3)
    kv4 = dt*Fv(v + kv3, o + ko3)
    ko4 = dt*Fo(o + ko3)

    nr = r + (kr1+2*kr2+2*kr3+kr4)/6
    nv = v + (kv1+2*kv2+2*kv3+kv4)/6
    no = o + (ko1+2*ko2+2*ko3+ko4)/6
    
    return nr, nv, no

for i in range(nsteps):
    rs[i], vs[i], os[i] = r, v, o
    r, v, o = SimRK4(r, v, o)

rT = rs.T
ax = plt.axes(projection='3d')
ax.plot3D(rT[0], rT[1], rT[2], 'gray')
plt.show()