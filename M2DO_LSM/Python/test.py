import numpy as np
import matplotlib.pyplot as plt
from py_lsmBind import py_LSM

nelx = 160
nely = 80

hole = np.array([[16, 14, 5],
                [48, 14, 5],
                [80, 14, 5],
                [112, 14, 5],
                [144, 14, 5],

                [32, 27, 5],
                [64, 27, 5],
                [96, 27, 5],
                [128, 27, 5],

                [16, 40, 5],
                [48, 40, 5],
                [80, 40, 5],
                [112, 40, 5],
                [144, 40, 5],

                [32, 53, 5],
                [64, 53, 5],
                [96, 53, 5],
                [128, 53, 5],

                [16, 66, 5],
                [48, 66, 5],
                [80, 66, 5],
                [112, 66, 5],
                [144, 66, 5]],dtype=np.float)

a = py_LSM(nelx,nely)
print(len(hole[:,0]))
a.add_holes(locx = hole[:,0], locy = hole[:,1], radius = hole[:,2])

a.set_levelset(True)

(bpts_xy, areafraction, seglength) = a.discretise()


plt.scatter(bpts_xy[:,0], bpts_xy[:,1], 10, seglength)
plt.show()


