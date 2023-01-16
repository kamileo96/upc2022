import numpy as np, matplotlib.pyplot as plt

dt = 0.01
maxt =10
g = 9.81

R = 0.11
m=0.4
lmb = 8*np.pi*1.6e-5*R*3/(m*2)
C_L = 0.3
C1 = C_L*R**3*np.pi*1.2
C_d = 0.25
C2 = C_d*1.2*np.pi*R**2/2
goalx = 3.66
goalz = 2.44
window = 0.75
nsteps = int(maxt/dt)
r0 = np.array([0., 0., R])

# rs = np.empty((nsteps, 3))
# vs = np.empty((nsteps, 3))
# os = np.empty((nsteps, 3))

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
    w=-3
    if r[1]<11 and r[2]<R: 
        w = -2
    if (abs(r[0])>=goalx-R or r[2]>=goalz-R) and r[1] >= 11:
        w=-1
    if (abs(goalx - abs(r[0])) <= R or abs(goalz - r[2]) <= R) and r[1] >= 11 and (abs(r[0])<goalx+R and r[2]<goalz+R):
        w=0
    if abs(r[0])<goalx-R and r[1] >= 11 and r[2]<goalz-R:
        w=1
    if w==1 and abs(r[0])>goalx-window+R and r[2]>goalz-window+R:
        w = 2
    if w==-3: print('Shot evaluation error')
    return w
def Shoot(r0, v0, o0):
    r, v, o = r0, v0, o0
    i=0
    while i<nsteps and r[2] >= R and r[1] <= 11:
        r, v, o = SimRK4(r, v, o)
        i+=1
        if(i>=nsteps): print("time out")
    return r

v_power = 35
o0 = np.array([0., 0., 0.])

precision = 4
numphis = 50*precision
LEFT_phi = 60#deg
RIGHT_phi = 120
UP_theta = 60
DOWN_theta = 90
phis = np.linspace(LEFT_phi/180*np.pi, RIGHT_phi/180*np.pi, numphis)[::-1]
numthetas = 30*precision
thetas = np.linspace(UP_theta/180*np.pi, DOWN_theta/180*np.pi, numthetas)
results = np.empty((numphis, numthetas))

for i in range(numphis):
    for j in range(numthetas):
        phi = phis[i]
        theta = thetas[j]
        v0 = v_power*np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        results[i][j] = ShotEvaluator(Shoot(r0, v0, o0))

im = plt.matshow(results.T)
im.set_extent(np.array([LEFT_phi, RIGHT_phi, DOWN_theta, UP_theta]))
boundaries=[-2.5, -1.5,-0.5,0.5,1.5,2.5]
ticks=[-2, -1, 0, 1, 2]
aspect=5
plt.colorbar(im, shrink=0.3, aspect=aspect, boundaries=boundaries, ticks=ticks)
plt.xlabel(r'$\theta $ [deg]')
plt.ylabel(r'$\varphi $ [deg]')
plt.show()

