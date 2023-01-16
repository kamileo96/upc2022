import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.colors import Normalize
dt = 0.001
maxt = 3
g = 9.81
distance = 11
m=0.4
R = 0.11
C_L = 0.3
C1 = C_L*R**3*np.pi*1.2
C_d = 0.25
C2 = C_d*1.2*np.pi*R**2/2
lmb = 8*np.pi*1.6e-5*R*3/(m*2)

goalx = 3.66
goalz = 2.44
window = 0.75
nsteps = int(maxt/dt)
r0 = np.array([0., 0., R])
r=r0
theta=np.pi*0.4
phi = np.pi*0.5
v = 18*np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
o = np.array([0., 0., 150.])

# v = np.array([7.2, 10., 7.9])
# o = np.array([0., 0., 10.])


# v = np.array([6., 10., 6.5])
# o = np.array([0., 0., 10.])

rs = np.empty((nsteps, 3))
vs = np.empty((nsteps, 3))
os = np.empty((nsteps, 3))

def Fdrag(v):
    return -C2*np.linalg.norm(v)*v

def Fmagnus(v, o):
    cr7 = np.cross(o, v)
    SIUU = np.linalg.norm(cr7)
    if SIUU == 0: return 0
    return C1*np.linalg.norm(v)*np.linalg.norm(o)*cr7/SIUU#/C_L*(0.8*R*np.linalg.norm(o)/np.linalg.norm(v)+0.12)

def fr(v):
    return v
def fv(v, o):
    fg = np.array([0., 0., -g])

    f = fg + Fdrag(v)/m + Fmagnus(v, o)/m
    return f
def fo(o):
    return -lmb*o

def SimRK4(r, v, o):
    kr1 = dt*fr(v)
    kv1 = dt*fv(v, o)
    ko1 = dt*fo(o)

    kr2 = dt*fr(v + kv1/2)
    kv2 = dt*fv(v + kv1/2, o + ko1/2)
    ko2 = dt*fo(o + ko1/2)

    kr3 = dt*fr(v + kv2/2)
    kv3 = dt*fv(v + kv2/2, o + ko2/2)
    ko3 = dt*fo(o + ko2/2)

    kr4 = dt*fr(v + kv3)
    kv4 = dt*fv(v + kv3, o + ko3)
    ko4 = dt*fo(o + ko3)

    nr = r + (kr1+2*kr2+2*kr3+kr4)/6
    nv = v + (kv1+2*kv2+2*kv3+kv4)/6
    no = o + (ko1+2*ko2+2*ko3+ko4)/6

    return nr, nv, no
def ShotEvaluator(r):
    w=-2
    if r[1]<distance and r[2]<R: 
        w = -1
    if (abs(r[0])>=goalx-R or r[2]>=goalz-R) and r[1] >= distance:
        w=0
    if abs(r[0])<goalx-R and r[1] >= distance and r[2]<goalz-R:
        w=1
    if w==1 and abs(r[0])>goalx-window+R and r[2]>goalz-window+R:
        w = 2
    if w==-2: print('Shot evaluation error')
    return w

goalieZ = goalz/2
def MinimalGoalkeeperVelocity(r, t):
    t = t-0.2
    vy = (r[2] - goalieZ +g*t**2/2)/t
    vx = np.linalg.norm(np.array([r[0], distance-r[1]]))/t
    return np.linalg.norm(np.array([vy, vx]))
goalvs = np.empty(nsteps)
mini = int(0.4/dt)
i=0
while i<nsteps and r[2] >= 0.9*R and r[1] <= distance:
    rs[i], vs[i], os[i] = r, v, o
    if i>=mini:
        goalvs[i] = MinimalGoalkeeperVelocity(r, dt*(i))
    r, v, o = SimRK4(r, v, o)
    i+=1

goalvs=goalvs[mini:i-1]
rs, vs, os = rs[:i-1], vs[:i-1], os[:i-1]
rT = rs.T
ax = plt.axes(projection='3d')
ax.plot3D(rT[0][:mini], rT[1][:mini], rT[2][:mini], 'blue')
normalize = Normalize(vmin=np.amin(goalvs), vmax=np.amax(goalvs))
goalvs=np.log(np.sqrt(goalvs))
goalvs = goalvs-np.amin(goalvs)
goalvs=goalvs/np.amax(goalvs)
ax.scatter(rT[0][mini:], rT[1][mini:], rT[2][mini:], c = plt.cm.viridis(goalvs)) 
GOAL = np.array([[-3.66, distance, 0.], [-3.66, distance, 2.44], [3.66, distance, 2.44], [3.66, distance, 0.]])
GOAL = GOAL.T
ax.plot3D(GOAL[0], GOAL[1], GOAL[2], 'gray')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
BOX = 12
ax.set_xlim3d([-BOX/2, BOX/2])
ax.set_ylim3d([0, BOX])
ax.set_zlim3d([0, BOX/2])

p = Rectangle((-goalx, goalz-window), window, window, color='green', alpha=0.5)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=distance, zdir="y")

p = Rectangle((goalx-window, goalz-window), window, window, color='green', alpha=0.5)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=distance, zdir="y")

# p = Rectangle((r0[0]-1.5, r0[2]), 3, 2, color='black', alpha=0.5)
# ax.add_patch(p)
# art3d.pathpatch_2d_to_3d(p, z=8.5, zdir="y")

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)*R+r0[0]
y = np.sin(u)*np.sin(v)*R+r0[1]
z = np.cos(v)*R+R
ax.plot_wireframe(x, y, z, color="r")

x = np.cos(u)*np.sin(v)*R+r[0]
y = np.sin(u)*np.sin(v)*R+r[1]
z = np.cos(v)*R+r[2]
ax.plot_wireframe(x, y, z, color="r")

ax.set_box_aspect(aspect = (1,1,0.5))
print(r)
print(ShotEvaluator(r))

plt.show()
print(r[1])
print(goalvs)
print(np.amin(goalvs), np.amax(goalvs))