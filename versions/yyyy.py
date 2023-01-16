import numpy as np, matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d


ax = plt.axes(projection='3d')
p = Rectangle((-1, 0), 2, 1)#, alpha =0.5)

ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=2, zdir="y")

u, v = np.mgrid[-1:1:20j, 0:1:10j]
y = u*0+2
ax.plot_wireframe(u, y, v, alpha=0.5, color='black')

#ax.quiver(0, np.sqrt(2)/2-0.5, np.sqrt(2)/2, 0, 1, 0, length=0.3, normalize=True, color='orange')
u, v = np.mgrid[0:np.pi:40j, 0:np.pi/2:20j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)-0.5
z = np.cos(v)
ax.plot_wireframe(x, y, z, alpha=0.5, color='black')
ax.set_xlim3d([-1.5, 1.5])
ax.set_ylim3d([-1, 2])
ax.set_zlim3d([0, 3])

ax.set_box_aspect(aspect = (1,1,1))

plt.show()