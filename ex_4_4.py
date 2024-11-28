"""
Ví dụ 4.4
Tính dầm đơn giản
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

I = 2510e-8; E = 210e11
ep = [E, I]

edofs=np.array([[1,2,3,4],     
                [3,4,5,6]])

ex = np.array([[0,3],
               [3,9]])

# ----------------------------- KHỐI 3 -----------------------------

Ke1=cfc.beam1e(ex[0],ep)
Ke2=cfc.beam1e(ex[1],ep)

K = np.zeros((6,6))          
f = np.zeros((6,1))

cfc.assem(edofs[0],K,Ke1)
cfc.assem(edofs[1],K,Ke2)

# ----------------------------- KHỐI 4 -----------------------------

bc=np.array([1,5])

f[2,0] = -10000

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed = cfc.extract_eldisp(edofs,d)

es1= cfc.beam1s(ex[0],ep,ed[0])
es2= cfc.beam1s(ex[1],ep,ed[1])

# ----------------------------- KHỐI 7 -----------------------------

ndis = d.reshape(3,2)
cfu.disp_array(ndis, ['uy','fi'])

cfu.disp_array(es1, ['Q_1','M_1'])
cfu.disp_array(es2, ['Q_2','M_2'])

# ------------------------------------------------------------------