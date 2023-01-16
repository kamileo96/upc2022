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
v0 = np.array([8., 10., 9.5])
o0 = np.array([0., 0., 10.])

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


v_power = 25
o0 = np.array([0., 0., 0.])

precision = 1
numphis = 50*precision
LEFT_phi = 60#deg
RIGHT_phi = 120
UP_theta = 60
DOWN_theta = 90
phis = np.linspace(LEFT_phi/180*np.pi, RIGHT_phi/180*np.pi, numphis)[::-1]
numthetas = 60
thetas = np.linspace(UP_theta/180*np.pi, DOWN_theta/180*np.pi, numthetas)
results = np.empty((numphis, numthetas))

correction = np.empty((numphis, numthetas, 3))

def corvector(r):
    r2 = r
    r2[1]=0
    #v1 = np.array([0., 0., goalz/2])-r2
    v1 = np.array([goalx-window/2, 0., goalz-window/2])-r2
    #u = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
    u = np.array([0., 1., 0.])
    v2 = np.cross(u, v1)
    v2 = v2/np.linalg.norm(v2)

    return v2*(1-np.exp(-np.linalg.norm(v1)**2/window**2))


for i in range(numphis):
    for j in range(numthetas):
        phi = phis[i]
        theta = thetas[j]
        v0 = v_power*np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        r = Shoot(r0, v0, o0)
        results[i][j] = ShotEvaluator(r)
        #if(results[i][j] < 2):
        #if(results[i][j] < 1):
        correction[i][j] = corvector(r)
        #else:
            #correction[i][j] = np.array([0.,0.,0.])
        

im = plt.matshow(results.T)
im.set_extent(np.array([LEFT_phi, RIGHT_phi, DOWN_theta, UP_theta]))
boundaries=[-2.5, -1.5,-0.5,0.5,1.5,2.5]
ticks=[-2, -1, 0, 1, 2]
aspect=5
plt.colorbar(im, shrink=0.3, aspect=aspect, boundaries=boundaries, ticks=ticks)
plt.show()


v_power = 25
o0 = np.array([0., 0., 0.])

r0 = np.array([0., 0., R])
numphis = 50*precision
LEFT_phi = 60#deg
RIGHT_phi = 120
UP_theta = 60
DOWN_theta = 90
phis = np.linspace(LEFT_phi/180*np.pi, RIGHT_phi/180*np.pi, numphis)[::-1]
numthetas = 60
thetas = np.linspace(UP_theta/180*np.pi, DOWN_theta/180*np.pi, numthetas)
results = np.empty((numphis, numthetas))
alphas = np.zeros((numphis, numthetas))
OMEGA = 70
for i in range(numphis):
    for j in range(numthetas):
        o0 = OMEGA*correction[i][j]
        phi = phis[i]
        theta = thetas[j]
        v0 = v_power*np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        r = Shoot(r0, v0, o0)
        w = ShotEvaluator(r)
        if w==-2:
            w=OMEGA/2
        elif w==-1:
            w=0
        elif w==0:
            w=OMEGA
        elif w>0:
            w = np.linalg.norm(o0)
            alphas[i][j] = 1
            
        results[i][j] = w
fig = plt.figure(frameon=False)
print(results.T.shape, alphas.shape)
im = plt.imshow(results.T)
im.set_cmap('hot_r')
im.set_extent(np.array([LEFT_phi, RIGHT_phi, DOWN_theta, UP_theta]))

boundaries=[-2.5, -1.5,-0.5,0.5,OMEGA]
ticks=[-2, -1, 0, 1, OMEGA/2]
aspect=5
#plt.colorbar(im, shrink=0.3, aspect=aspect, boundaries=boundaries, ticks=ticks)

im2 = plt.imshow(results.T)
im2.set_cmap('nipy_spectral')
im2.set_alpha(alphas.T)

im2.set_extent(np.array([LEFT_phi, RIGHT_phi, DOWN_theta, UP_theta]))
plt.show()
