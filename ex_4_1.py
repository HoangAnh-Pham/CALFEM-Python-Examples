"""
Ví dụ 4.1
Tính thanh chịu kéo có tiết diện thay đổi
"""
# ----------------------------- KHỐI 1 -----------------------------

import calfem.core as cfc
import calfem.utils as cfu
import numpy as np

# ----------------------------- KHỐI 2 -----------------------------

A1,E1,A2,E2,L,P = 6e-4, 2e8, 3e-4, 2e8, 1.6, 50
ep1, ep2 = [E1,A1], [E2,A2]
ex1, ex2 = [0,L], [L,2*L]
edofs = np.array([[1,2],[2,3]])

# ----------------------------- KHỐI 3 -----------------------------

Ke1 = cfc.bar1e(ex1,ep1)
Ke2 = cfc.bar1e(ex2,ep2)

K = np.zeros((3,3))
f = np.zeros((3,1))
cfc.assem(edofs[0],K,Ke1)
cfc.assem(edofs[1],K,Ke2)

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1])
f[2] = P

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed = cfc.extract_eldisp(edofs,d)
es1 = cfc.bar1s(ex1,ep1,ed[0])
es2 = cfc.bar1s(ex2,ep2,ed[1])

# ----------------------------- KHỐI 7 -----------------------------


cfu.disp_array(d,['Chuyển vị nút'])
cfu.disp_array(r,['Phản lực nút'])
cfu.disp_array(es1,['N1'])
cfu.disp_array(es2,['N2'])

# ------------------------------------------------------------------