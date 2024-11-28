"""
Ví dụ 4.8
Tính hệ liên hợp
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----------------------------- KHỐI 2 -----------------------------

A1, I1, A2 = 4e-3, 5.4e-5, 1e-3  
E = 2e11
q = 10

ep1 = [E,A1,I1]
ep2 = [E,A2]

edof1 = np.array([[1, 2, 3, 4, 5, 6], 
                  [4, 5, 6, 7, 8, 9],
                  [7, 8, 9, 10, 11, 12]])

edof2 = np.array([[13, 14, 4, 5],     
                  [13, 14, 7, 8]])

L, H = 2, 2
coords = np.array([[0, H],
                   [L, H],
                   [2*L, H],
                   [3*L, H],
                   [0, 0]])

dofs = np.array([[ 1,  2,  3],
                 [ 4,  5,  6],
                 [ 7,  8,  9],
                 [10, 11, 12],
                 [13, 14, 15]])

eq = np.array([[0,0],
               [0,-q],
               [0,-q]])

ex1, ey1 = cfc.coordxtr(edof1, coords, dofs, 2)
ex2, ey2 = cfc.coordxtr(edof2, coords, dofs, 2)

nod = np.max(dofs)
noe1 = len(edof1)
noe2 = len(edof2)

# ----------------------------- KHỐI 3 -----------------------------

K = np.zeros((nod,nod))      
f = np.zeros((nod,1)) 

for i in range(noe1):
    Ke, fe = cfc.beam2e(ex1[i],ey1[i],ep1,eq[i])    
    cfc.assem(edof1[i],K,Ke,f,fe)
for i in range(noe2):
    Ke = cfc.bar2e(ex2[i],ey2[i],ep2)    
    cfc.assem(edof2[i],K,Ke)

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1,2,3,13,14,15])

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed1 = cfc.extract_eldisp(edof1,d)
ed2 = cfc.extract_eldisp(edof2,d)

nsec = 5
es1 = np.zeros((noe1,nsec,3))
es2 = np.zeros((noe2,2,1))

for i in range(noe1):
    es1[i],edi,eci = cfc.beam2s(
        ex1[i],ey1[i],ep1,ed1[i],eq[i],nsec)
for i in range(noe2):
    es2[i] = cfc.bar2s(ex2[i],ey2[i],ep2,ed2[i])

# ----------------------------- KHỐI 7 -----------------------------

cfu.disp_array(es1[0],['N_1','Q_1','M_1'])
cfu.disp_array(es1[1],['N_2','Q_2','M_2'])
cfu.disp_array(es1[2],['N_3','Q_3','M_3'])

cfu.disp_array(es2[0],['N_4'])
cfu.disp_array(es2[1],['N_5'])

# ------------------------------------------------------------------