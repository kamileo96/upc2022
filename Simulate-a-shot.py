import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d

#constants in SI units
dt = 0.01 #can be 0.001
maxt = 3
g = 9.81
m=0.4
R = 0.11
C_L = 0.3 #Lift coefficient 
C1 = C_L*R**3*np.pi*1.2 #effective lift coeff.
C_d = 0.25 #drag coeff.
C2 = C_d*1.2*np.pi*R**2/2 #eff. drag
lmb = 8*np.pi*1.6e-5*R*3/(m*2) #lambda for omega dissipation, 
#uses Reynolds number

goalx = 3.66 #goal size / 2
goalz = 2.44
window = 0.75 #window square size
nsteps = int(maxt/dt)


#################################################
#EASY CHANGE:
distance = 11 #of the goal. can be changed i. e. 
#for free kicks (ball starts at 0. 0.)

xoffset = 0
#initial velocity:
v0 = 20 #m/s
theta=90/180*np.pi #up to 180 deg
phi = 90/180*np.pi

o = np.array([0., 0., 70.]) #\omega in rad/s
################################################



v = v0*np.array([np.sin(theta)*np.cos(phi), 
np.sin(theta)*np.sin(phi), np.cos(theta)])

r0 = np.array([0., xoffset, R])
r=r0 
rs = np.empty((nsteps, 3))
vs = np.empty((nsteps, 3))
os = np.empty((nsteps, 3))

def Fdrag(v):
    return -C2*np.linalg.norm(v)*v

def Fmagnus(v, o):
    cr7 = np.cross(o, v)
    SIUU = np.linalg.norm(cr7)
    if SIUU == 0: return 0
    return C1*np.linalg.norm(v)*np.linalg.norm(o)*cr7/SIUU

def fr(v):
    return v

def fv(v, o):
    fg = np.array([0., 0., -g])

    f = fg + Fdrag(v)/m + Fmagnus(v, o)/m
    return f
def fo(o):
    return -lmb*o

def SimRK4(r, v, o):
    #Runge Kuttas method of the 4th order
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
    #-2 error
    #-1 hit the ground
    # 0 hit the frame 
    # 1 hit the goal but not
    # 2 hit the upper corner!!
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


"""
goalieZ = goalz/2
def MinimalGoalkeeperVelocity(r, t):
    t = t-0.2 #goalie reaction time
    vy = (r[2] - goalieZ +g*t**2/2)/t
    vx = np.linalg.norm(np.array([r[0], distance-r[1]]))/t
    return np.linalg.norm(np.array([vy, vx]))
goalvs = np.empty(nsteps)
mini = int(0.25/dt) #when the scatter starts
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
ax.scatter(rT[0][mini:], rT[1][mini:], rT[2][mini:], 
c = plt.cm.viridis(goalvs/np.amax(goalvs))) """
#^use this to emulate a goalkeeper


########
i=0
while i<nsteps and r[2] >= 0.99*R and r[1] <= distance:
    rs[i], vs[i], os[i] = r, v, o
    r, v, o = SimRK4(r, v, o)
    i+=1


rs, vs, os = rs[:i-1], vs[:i-1], os[:i-1]
rT = rs.T
ax = plt.axes(projection='3d')
ax.plot3D(rT[0], rT[1], rT[2], 'blue')

##########
GOAL = np.array([[-goalx, distance, 0.], [-goalx, distance, goalz], 
[goalx, distance, goalz], [goalx, distance, 0.]])
GOAL = GOAL.T
ax.plot3D(GOAL[0], GOAL[1], GOAL[2], 'gray')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
########### Here you can expand the viewing box, but adjust the aspect!
BOX = 12
ax.set_xlim3d([-BOX/2, BOX/2])
ax.set_ylim3d([0, BOX])
ax.set_zlim3d([0, BOX/2])
ax.set_box_aspect(aspect = (1,1,0.5))


p = Rectangle((-goalx, goalz-window), 
window, window, color='green', alpha=0.5)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=distance, zdir="y")

p = Rectangle((goalx-window, goalz-window), 
window, window, color='green', alpha=0.5)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=distance, zdir="y")


##### Here you can add a wall
# p = Rectangle((r0[0]-1.5, r0[2]), 3, 2, 
# color='black', alpha=0.5)
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

plt.show()