"""
Ví dụ 4.7
Tính hệ khung phẳng với thanh liên kết khớp
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

b1, h1, b2, h2= 0.2, 0.4, 0.2, 0.2
E = 2e7
q, P, M = 10, 50, 100

A1 = b1*h1; I1 = b1*h1**3/12
A2 = b2*h2; I2 = b2*h2**3/12
ep = np.array([[E,A1,I1],
               [E,A1,I1],
               [E,A1,I1],
               [E,A2,I2],
               [E,A2,I2],
               [E,A2,I2]])

edofs=np.array([[ 7,  8,  9, 19, 20, 21],
                [19, 20, 21, 10, 11, 12],
                [10, 11, 12, 13, 14, 15],   
                [ 1,  2,  3, 16, 17, 18],
                [16, 17, 18,  7,  8,  9],
                [ 4,  5,  6, 10, 11, 12]])

L, H = 3, 2
coords = np.array([[0, 0],
                   [2*L, 0],
                   [0, 2*H],
                   [2*L, 2*H],
                   [4*L, 2*H],
                   [0, H],
                   [L, 2*H]])

dofs = np.array([[ 1,  2,  3],
                 [ 4,  5,  6],
                 [ 7,  8,  9],
                 [10, 11, 12],
                 [13, 14, 15],
                 [16, 17, 18],
                 [19, 20, 21]])

eq = np.zeros((len(edofs),2))
eq[2] = np.array([0,-q])

ex, ey = cfc.coordxtr(edofs, coords, dofs, 2)

edofs[2,2] = 22

nod = np.max(edofs)
noe = len(edofs)
non = len(dofs)
# ----------------------------- KHỐI 3 -----------------------------

K = np.zeros((nod,nod))      
f = np.zeros((nod,1)) 

for i in range(noe):
    Ke, fe = cfc.beam2e(ex[i],ey[i],ep[i],eq[i])    
    cfc.assem(edofs[i],K,Ke,f,fe)

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1,2,3,4,5,6,13,14,15])

f[15,0] = f[15,0] + P
f[8,0] = f[8,0] - M
f[19,0] = f[19,0] - P
f[10,0] = f[10,0] - P

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

eds = cfc.extract_eldisp(edofs,d)

nsec = 11
es = np.zeros((noe,nsec,3))
ed = np.zeros((noe,nsec,2)) 

for i in range(noe):
    es[i],ed[i],eci = cfc.beam2s(
        ex[i],ey[i],ep[i],eds[i],eq[i],nsec)

# ----------------------------- KHỐI 7 -----------------------------

cfu.disp_array(es[2],['N_3','Q_3','M_3'])

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