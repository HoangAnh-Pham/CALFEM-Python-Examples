"""
Ví dụ 4.6
Tính hệ khung phẳng
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

b1, h1, b2, h2= 0.2, 0.2, 0.2, 0.3
E = 2e7
q, P, M = 10, 50, 100

A1 = b1*h1; I1 = b1*h1**3/12
A2 = b2*h2; I2 = b2*h2**3/12
ep = np.array([[E,A1,I1],
               [E,A2,I2],
               [E,A2,I2]])

edofs = np.array([[ 1,  2,  3,  4,  5,  6],
                  [ 7,  8,  9, 10, 11, 12],
                  [10, 11, 12,  1,  2,  3]])

L, H = 3, 2
coords = np.array([[L, 2*H],
                   [2*L, 2*H],
                   [0, 0],
                   [L/2, H]])

dofs = np.array([[ 1,  2,  3],
                 [ 4,  5,  6],
                 [ 7,  8,  9],
                 [10, 11, 12]])

eq = np.array([[0,-q],       
               [0,0],
               [0,0]])

ex, ey = cfc.coordxtr(edofs, coords, dofs, 2)

nod = np.max(dofs)
noe = len(edofs)
non = len(dofs)
# ----------------------------- KHỐI 3 -----------------------------

K = np.zeros((nod,nod))      
f = np.zeros((nod,1)) 

for i in range(noe):
    Ke, fe = cfc.beam2e(ex[i],ey[i],ep[i],eq[i])    
    cfc.assem(edofs[i],K,Ke,f,fe)

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([4,5,6,7,8])

f[10,0] = f[10,0] - P
f[2,0] = f[2,0] - M

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

eds = cfc.extract_eldisp(edofs,d)

nsec = 5
es = np.zeros((noe,nsec,3))
ed = np.zeros((noe,nsec,2)) 

for i in range(noe):
    es[i],ed[i],eci = cfc.beam2s(
        ex[i],ey[i],ep[i],eds[i],eq[i],nsec)

# ----------------------------- KHỐI 7 -----------------------------

cfu.disp_array(es[0],['N_1','Q_1','M_1'])

cfv.figure(1)
cfv.eldraw2(ex, ey)
for exi, eyi, edi in zip(ex, ey, ed):
    cfv.dispbeam2(
        exi, eyi, edi, plotpar=[1,4,1], sfac=1e2)
cfv.title('Sơ đồ biến dạng')
cfv.axis('equal')

plotpar = [4, 2]
sfac = 0.05

cfv.figure(2)
for exi, eyi, esi in zip(ex, ey, es):
    cfv.secforce2(exi, eyi, -esi[:,0], plotpar, sfac)
cfv.title('Biểu đồ lực dọc')
cfv.axis('equal'); 

cfv.figure(3)
for exi, eyi, esi in zip(ex, ey, es):
    cfv.secforce2(exi, eyi, esi[:,1], plotpar, sfac)
cfv.title('Biểu đồ lực cắt')
cfv.axis('equal'); 

cfv.figure(4)
for exi, eyi, esi in zip(ex, ey, es):
    cfv.secforce2(exi, eyi, esi[:,2], plotpar, sfac)
cfv.title('Biểu đồ mômen uốn')
cfv.axis('equal'); 

# ------------------------------------------------------------------