'''
Ví dụ 3.1: 
Tính hệ 3 lò xo
'''
# ----------------------------- KHỐI 1 -----------------------------

# Nạp các thư viện
import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----------------------------- KHỐI 2 -----------------------------

# Lập ma trận kết nối phần tử  
edofs = np.array([
    [1, 2],      # lò xo 1 nối chuyển vị 1 và 2
    [2, 3],      # lò xo 2 nối chuyển vị 2 và 3
    [2, 3],      # lò xo 3 nối chuyển vị 2 và 3
])

# ----------------------------- KHỐI 3 -----------------------------

# Tính ma trận độ cứng của lò xo
k = 1500.
Ke1 = cfc.spring1e(k)
Ke2 = cfc.spring1e(2*k)   
# Tạo ma trận độ cứng và véc tơ lực nút của hệ
K = np.zeros((3,3))
f = np.zeros((3,1))
# Ghép nối các lò xo
cfc.assem(edofs[0], K, Ke2)
cfc.assem(edofs[1], K, Ke1)
cfc.assem(edofs[2], K, Ke2)

# ----------------------------- KHỐI 4 -----------------------------

# Khai báo chuyển vị bị ngăn cản
bc = np.array([1, 3])
# Gán tải trọng nút tương ứng chuyển vị 2
f[1] = 100.0

# ----------------------------- KHỐI 5 -----------------------------

# Giải hệ phương trình cân bằng toàn hệ 
d, r = cfc.solveq(K, f, bc)
# In các kết quả tính toán
cfu.disp_array(d,['Chuyển vị nút'])
cfu.disp_array(r,['Phản lực nút'])

# ----------------------------- KHỐI 6 -----------------------------

# Trích xuất chuyển vị nút của phần tử
ed = cfc.extract_ed(edofs, d)
# Tính lực dọc trong các lò xo
es1 = cfc.spring1s(2*k, ed[0])
es2 = cfc.spring1s(k, ed[1])
es3 = cfc.spring1s(2*k, ed[2])

# ----------------------------- KHỐI 7 -----------------------------

# Hiển thị lực dọc trong lò xo
print('N1 = '+str(es1))
print('N2 = '+str(es2))
print('N3 = '+str(es3))

# ------------------------------------------------------------------