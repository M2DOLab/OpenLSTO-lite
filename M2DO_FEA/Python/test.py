from pyBind import py_FEA, py_Sensitivity
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg

from matching_dofIds import matching_dof, matching_el

nelx = 10
nely = 5

a = py_FEA(nelx = nelx, nely = nely, element_order = 2)
[node, elem, elem_dof] = a.get_mesh()
a.set_material(1.0,0.3,1)

# plt.figure(1)
# plt.axis([0,max(node[:,0]),0,max(node[:,1])])
# elem_ext = np.zeros([elem.shape[0],elem.shape[1]+1],dtype = int)
# elem_ext[:,0:4] = elem
# elem_ext[:,4] = elem[:,0]

coord = np.array([0,0])
tol = np.array([1e-3,1e10])

a.set_boundary(coord = coord,tol = tol)
BCid = a.get_boundary()

coord = np.array([nelx,nely/2])
tol = np.array([1,1]) 
la = a.set_force(coord = coord,tol = tol, direction = 1, f = -1) 
la = np.array(la)
(rows, cols, vals) = a.compute_K()

u = a.solve_FE()

print min(u)

nprows = np.array(rows)
npcols = np.array(cols)
npvals = np.array(vals)

nELEM = elem.shape[0]
nNODE = node.shape[0]
nDOF  = nNODE * 2
nDOF_withLag  = nNODE * 2 + len(BCid)

mtx = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF_withLag, nDOF_withLag))

# force in np
# rhs = np.zeros((nNODE,2))

# fixedNodes = (abs(node[:,0]-coord[0]) < tol[0]) & (abs(node[:,1]-coord[1]) < tol[1])
# nid = np.where(fixedNodes == True)[0]

# rhs[nid,1] = -1
GF = np.zeros(nDOF_withLag)
GF[:nDOF] = la #rhs.transpose().flatten(order = 'F')


sol = scipy.sparse.linalg.spsolve(mtx, GF) # direct solver is used

print min(sol)

# order of dof is not conventional...

# b = py_Sensitivity(a,state)
# gpts_sens = b.compute_compliance_sens()

# bpts_sens = b.comptue_boundary_sens(np.array([[10,5],[10,10]]))

# (rowss, colss, valss) = a.compute_K_SIMP(np.ones(nELEM))
# nprows = np.array(rowss)
# npcols = np.array(colss)
# npvals = np.array(valss)
# mtxs = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF, nDOF))
# states = scipy.sparse.linalg.spsolve(mtxs, GF) # direct solver is used


# (rows1, cols1, vals1) = a.compute_K_PARM(np.ones(nNODE))
# nprows = np.array(rows1)
# npcols = np.array(cols1)
# npvals = np.array(vals1)

# mtx1 = scipy.sparse.csc_matrix((npvals, (nprows, npcols)), shape = (nDOF, nDOF))
# state1 = scipy.sparse.linalg.spsolve(mtx1, GF) # direct solver is used

# plt.plot(state,'bo')
# plt.plot(states,'r',)
# plt.plot(state1,'k')
# plt.axis('equal'), plt.legend(['0','SIMP','PARM'])
# plt.grid('on'), plt.show()

# (rows_ds, cols_ds, vals_ds) = a.compute_K_SIMP_derivs(state)
# # (rows_d1, cols_d1, vals_d1) = a.compute_K_PARM_derivs(state)

# nprows_ds = np.array(rows_ds)
# nprows_ds = np.array(rows_ds)
# npcols_ds = np.array(cols_ds)

# ## matching dofs
# (m2doIDs, m2do2fem2d) = matching_dof()
# (e_m2doIDs, e_m2do2fem2d) = matching_el()
# rows_ds_new = np.zeros(nprows_ds.shape, dtype=int)
# cols_ds_new = np.zeros(npcols_ds.shape, dtype=int)

# # rows_new = np.zeros(nprows.shape, dtype=int)
# states_new = np.zeros(states.shape, dtype=float)

# for mm in range(len(rows_ds)):
#     id0 = np.where(m2doIDs==rows_ds[mm])[0]
#     id1 = m2do2fem2d[id0]
#     rows_ds_new[mm] = id1

# for mm in range(len(cols_ds)):
#     id0 = np.where(e_m2doIDs==cols_ds[mm])[0]
#     id1 = e_m2do2fem2d[id0]
#     cols_ds_new[mm] = id1


#     # id0 = np.where(m2doIDs==cols_ds[mm])[0]
#     # id1 = m2do2fem2d[id0]
#     # cols_ds_new[mm] = id1
#     # print((mm,rows_ds[mm],cols_ds[mm],rows_ds_new[mm],cols_ds_new[mm]))
# # states_new[m2do2fem2d] = states[m2doIDs]

# # cols_ds_new = npcols_ds
# npvals_ds = np.array(vals_ds)
# # print(npvals_ds)

# # sortid = cols_ds_new.argsort().astype(int)
# # print(cols_ds_new)
# # rows_ds_new = rows_ds_new[sortid]
# # npvals_ds = npvals_ds[sortid]
# # plt.scatter(rows_ds_new, cols_ds_new, 10, npvals_ds)
# # plt.scatter(np.array(rows_ds), np.array(cols_ds), 10, np.array(npvals_ds))

# # print(rows_ds_new.shape)
# # print(cols_ds_new.shape)
# # print(npvals_ds.shape)
# # print(nDOF)
# # print(np.array(cols_ds))    
# # print(cols_ds_new)

# mtxs = scipy.sparse.csc_matrix((npvals_ds, (rows_ds_new, cols_ds_new)), shape = (nDOF, nELEM))
# mtx_dens = mtxs.todense()
# mtx2 = scipy.sparse.csc_matrix(mtx_dens)

# aa = np.zeros(nDOF)
# for ii in range(nDOF):
#     aa[ii] = np.sum(mtx_dens[ii,:])

# plt.plot(aa)
# plt.axis("tight")
# plt.show()

# print(aa)
# # plt.subplot(2,1,1)
# # plt.plot(rows_ds_new, npvals_ds,'o')
# # plt.subplot(2,1,2)
# # plt.plot(cols_ds_new, npvals_ds,'o')
# # # plt.plot(nprows_ds,rows_ds_new)
# # plt.show()
