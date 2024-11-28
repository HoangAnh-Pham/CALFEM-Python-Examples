"""
Ví dụ 4.2
Tính dàn phẳng 2 thanh
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

S1 = 0.2*0.2; S2 = 0.2*0.3; E = 2.e7
ep1, ep2 = [E,S1], [E,S2]

edofs = np.array([[3,4,5,6],
                  [1,2,5,6]])
coords = np.array([[0, 0],
                   [0, 4.],
                   [3., 4.]])
dofs = np.array([[1, 2],
                 [3, 4],
                 [5, 6]])
ex, ey = cfc.coordxtr(edofs, coords, dofs, 2)

# ----------------------------- KHỐI 3 -----------------------------

Ke1 = cfc.bar2e(ex[0],ey[0],ep1)
Ke2 = cfc.bar2e(ex[1],ey[1],ep2)

K = np.zeros((6,6))
f = np.zeros((6,1))
cfc.assem(edofs[0],K,Ke1)
cfc.assem(edofs[1],K,Ke2)

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1,2,3,4])
f[5] = -50

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed = cfc.extract_eldisp(edofs,d)
es1 = cfc.bar2s(ex[0],ey[0],ep1,ed[0])
es2 = cfc.bar2s(ex[1],ey[1],ep2,ed[1])

# ----------------------------- KHỐI 7 -----------------------------

cfv.eldraw2(ex, ey)
cfv.title('Sơ đồ hình học')

cfu.disp_array(d[dofs[2]-1],['Chuyển vị 5 và 6'])
cfu.disp_array(es1,['N1'])
cfu.disp_array(es2,['N2'])


# ------------------------------------------------------------------