"""
Ví dụ 4.5
Tính dầm liên tục
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

b, h, E, k = 0.2, 0.3, 2e7, 1e4
q, P, M = 10, 50, 100
I = b*h**3/12; 
ep = [E, I]

edofs = np.array([[1, 2, 3, 4],
                  [3, 4, 5, 6],
                  [5, 6, 7, 8],
                  [7, 8, 9, 10]])

coords = np.array([[0, 0],
                   [2, 0],
                   [4, 0],
                   [9, 0],
                   [13, 0]])

dofs = np.array([[1, 2],
                 [3, 4],
                 [5, 6],
                 [7, 8],
                 [9, 10]])

ex, ey = cfc.coordxtr(edofs, coords, dofs, 2)

nod = np.max(dofs)
noe = len(edofs)
non = len(dofs)

eq = np.array([[0],[0],[0],[-q]])

# ----------------------------- KHỐI 3 -----------------------------

K = np.zeros((nod,nod))      
f = np.zeros((nod,1)) 

for i in range(noe):
    Ke, fe = cfc.beam1e(ex[i],ep,eq[i])    
    cfc.assem(edofs[i],K,Ke,f,fe)

K[4,4] = K[4,4] + k

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1,2,7,9])

f[2,0] = f[2,0] - P
f[4,0] = f[4,0] - P
f[7,0] = f[7,0] - M

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed = cfc.extract_eldisp(edofs,d)

nsec = 5
es = np.zeros((noe,nsec,2))
for i in range(noe):
    es[i],edi,eci = cfc.beam1s(
        ex[i],ep,ed[i],eq[i],nsec)

# ----------------------------- KHỐI 7 -----------------------------

plotpar = [4, 2]
sfac = 1

# Vẽ biểu đồ lực cắt
cfv.figure(1)
for exi, eyi, esi in zip(ex, ey, es):
    cfv.secforce2(exi, eyi, esi[:,0], plotpar, sfac)
cfv.title('Biểu đồ lực cắt')

# Vẽ biểu đồ lực cắt
cfv.figure(2)
for exi, eyi, esi in zip(ex, ey, es):
    cfv.secforce2(exi, eyi, esi[:,1], plotpar, sfac)
cfv.title('Biểu đồ mômen uốn')

# ------------------------------------------------------------------