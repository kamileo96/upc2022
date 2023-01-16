import numpy as np
C1 = 1

v = np.array([0.,1.,0.])
o = np.array([0.,0.,1.])
def Fmagnus(v, o):
    cr7 = np.cross(o, v)
    return C1*np.dot(v, v)*cr7/np.linalg.norm(cr7)

print(Fmagnus(v, o))
print(np.cross(np.array([2, 0, 0]), np.array([0,2,0])))