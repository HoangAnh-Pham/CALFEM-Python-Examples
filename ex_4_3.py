"""
Ví dụ 4.3
Tính hệ dàn phẳng 9 phần tử sử dụng vòng lặp
"""
# ----------------------------- KHỐI 1 -----------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----------------------------- KHỐI 2 -----------------------------

S1 = 0.2*0.2; E = 2e7; k = 10000
ep = [E,S1]

edofs=np.array([[ 1,  2,  3,  4],
                [ 3,  4,  5,  6],
                [ 1,  2,  7,  8],
                [ 7,  8, 11, 12],
                [ 7,  8,  3,  4],
                [ 3,  4, 11, 12],
                [ 3,  4,  9, 10],
                [ 9, 10,  5,  6],
                [11, 12,  9, 10]])
a, b = 3., 2.
coords = np.array([[0, 0],
                   [2*a, 0],
                   [4*a, 0],
                   [a, b],
                   [3*a, b],
                   [2*a, 2*b]])
dofs = np.array([[1, 2],
                 [3, 4],
                 [5, 6],
                 [7, 8],
                 [9, 10],
                 [11, 12]])
ex, ey = cfc.coordxtr(edofs, coords, dofs, 2)

nod = np.max(dofs)
noe = len(edofs)
non = len(dofs)

# ----------------------------- KHỐI 3 -----------------------------

K = np.zeros((nod,nod))
f = np.zeros((nod,1))

for i in range(noe):
    Ke = cfc.bar2e(ex[i],ey[i],ep)
    cfc.assem(edofs[i],K,Ke) 
    
# Thêm độ cứng gối đàn hồi k=10000 
K[3,3] = K[3,3] + k

# ----------------------------- KHỐI 4 -----------------------------

bc = np.array([1,2,6])

# Gán tải trọng tập trung -50 ứng với chuyển vị 8,10,12
f[7,0] = -50
f[9,0] = -50
f[11,0] = -50
# Gán tải trọng tập trung 50 ứng với chuyển vị 11
f[10,0] = 50

# ----------------------------- KHỐI 5 -----------------------------

d, r = cfc.solveq(K,f,bc)

# ----------------------------- KHỐI 6 -----------------------------

ed = cfc.extract_eldisp(edofs,d)
for i in range(noe):
    es = cfc.bar2s(ex[i],ey[i],ep,ed[i])
    print('Lực dọc',i+1,':',es[0])

# ----------------------------- KHỐI 7 -----------------------------

cfv.eldraw2(ex, ey)
cfv.title('Sơ đồ hình học')

ndis = d.reshape(non,2)
cfu.disp_array(ndis, ['ux','uy'])

# ------------------------------------------------------------------