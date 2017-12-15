from pyBind import py_FEA, py_Sensitivity
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg

nelx = 10
nely = 10

a = py_FEA(nelx = nelx, nely = nely, element_order = 2)
[node, elem] = a.get_mesh()
a.set_material(1.0,0.3,1)

plt.figure(1)
plt.axis([0,max(node[:,0]),0,max(node[:,1])])
elem_ext = np.zeros([elem.shape[0],elem.shape[1]+1],dtype = int)
elem_ext[:,0:4] = elem
elem_ext[:,4] = elem[:,0]

coord = np.array([0,0])
tol = np.array([1e-3,1e10])

a.set_boundary(coord = coord,tol = tol)


coord = np.array([nelx,nely/2])
tol = np.array([1,1])
 
a.set_force(coord = coord,tol = tol, direction = 1, f = -1) 

(rows, cols, vals) = a.compute_K()
nprows = np.array(rows)
npcols = np.array(cols)
npvals = np.array(vals)

nELEM = elem.shape[0]
nNODE = node.shape[0]
nDOF  = nNODE * 2

# print ([rows, cols, vals])
mtx = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF, nDOF))

def mv(v):
    return array([ 2*v[0], 3*v[1]])

# force in np
rhs = np.zeros((nNODE,2))

fixedNodes = (abs(node[:,0]-coord[0]) < tol[0]) & (abs(node[:,1]-coord[1]) < tol[1])
nid = np.where(fixedNodes == True)[0]

rhs[nid,1] = -1
GF = rhs.transpose().flatten(order = 'F')

state = scipy.sparse.linalg.spsolve(mtx, GF) # direct solver is used
# order of dof is not conventional...

b = py_Sensitivity(a,state)
gpts_sens = b.compute_compliance_sens()

bpts_sens = b.comptue_boundary_sens(np.array([[10,5],[10,10]]))

(rowss, colss, valss) = a.compute_K_SIMP(np.ones(nELEM))
nprows = np.array(rowss)
npcols = np.array(colss)
npvals = np.array(valss)
mtxs = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF, nDOF))
states = scipy.sparse.linalg.spsolve(mtxs, GF) # direct solver is used

(rows1, cols1, vals1) = a.compute_K_PARM(np.ones(nNODE))
nprows = np.array(rows1)
npcols = np.array(cols1)
npvals = np.array(vals1)

mtx1 = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF, nDOF))
state1 = scipy.sparse.linalg.spsolve(mtx1, GF) # direct solver is used

(rows_ds, cols_ds, vals_ds) = a.compute_K_SIMP_derivs(state)
(rows_d1, cols_d1, vals_d1) = a.compute_K_PARM_derivs(state)

plt.plot(state,'bo')
plt.plot(states,'r',)
plt.plot(state1,'k')
plt.axis('equal'), plt.legend(['0','SIMP','PARM'])
plt.grid('on'), plt.show()
