'''
Ví dụ 5.5
Mô hình vách có lỗ mở
'''
# ------------------------------------------------------------------
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import calfem.geometry as cfg
import calfem.mesh as cfm
import numpy as np

cfv.close_all()
# ------------------------------------------------------------------
# Tạo dạng hình học
g = cfg.Geometry()
bottom_edge = 20
left_edge = 21

A, H, a, h = 3.2, 3.4, 1.2, 2.4
g.point([0, 0])     # point 0
g.point([A, 0])     # point 1
g.point([0, H])     # point 2
g.point([A, H])     # point 3
g.point([0, 2*H])   # point 4
g.point([A, 2*H])   # point 5

g.point([1, 0])     # point 6     
g.point([1+a, 0])   # point 7 
g.point([1, h])     # point 8
g.point([1+a, h])   # point 9

g.spline([0, 6], marker=bottom_edge)    # line 0
g.spline([7, 1], marker=bottom_edge)    # line 1
g.spline([4, 5])                        # line 2
g.spline([0, 2], marker=left_edge)      # line 3
g.spline([1, 3])                        # line 4
g.spline([2, 4], marker=left_edge)      # line 5
g.spline([3, 5])                        # line 6

g.spline([8, 9])                        # line 7
g.spline([6, 8])                        # line 8
g.spline([7, 9])                        # line 9

g.point([1, H])     # point 10
g.point([1+a, H])   # point 11
g.point([1, H+h])   # point 12
g.point([1+a, H+h]) # point 13

g.spline([10, 11])                      # line 10
g.spline([12, 13])                      # line 11
g.spline([10, 12])                      # line 12
g.spline([11, 13])                      # line 13

g.surface([0, 8, 7, 9, 1, 4, 6, 2, 5, 3],
          [[10,13,11,12]])

cfv.draw_geometry(g); 
cfv.axis('off')

# ------------------------------------------------------------------
# Tự động chia lưới
mesh = cfm.GmshMesh(g)
mesh.el_type = 2            # Element type
mesh.dofs_per_node = 2      # Degrees of freedom per node.
mesh.el_size_factor = 0.2   # Factor that changes element sizes.

coords, edofs, dofs, bdofs, elementmarkers = mesh.create()

ex, ey = cfc.coordxtr(edofs, coords, dofs)

cfv.figure()
cfv.draw_mesh(
    coords=coords,
    edof=edofs,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title='Lưới phần tử'
)

# ------------------------------------------------------------------
# Số liệu tính toán
E,nu,t,q = 3e6, 0.25, 0.1, 50.
ptype = 1
ep = [ptype,t]
D = cfc.hooke(ptype,E,nu)

# lập ma trận độ cứng của hệ
nod = np.size(dofs)
K = np.zeros([nod, nod])
f = np.zeros([nod,1])

for eldof, elx, ely in zip(edofs, ex, ey):
    Ke = cfc.plante(elx, ely, ep, D)
    cfc.assem(eldof, K, Ke)

# Gán điều kiện biên ngàm ở cạnh đáy
bc = np.array([],'i')
bcVal = np.array([],'f')    
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 
                        bottom_edge, 0.0, 0)

# Gán tải trọng phân bố trên cạnh trái
lc = np.array([],'i')
lcVal = np.array([],'f')

ndof = mesh.dofs_per_node
no_lpoint = len(bdofs[left_edge])/ndof
q_point = q*2*H/no_lpoint
lc, lcVal = cfu.applybc(bdofs, lc, lcVal, 
                        left_edge, q_point, 1)
f[lc-1,0] = lcVal

# Giải hệ phương trình
d, r = cfc.solveq(K,f,bc)

# Tính toán ứng suất
ed = cfc.extract_eldisp(edofs,d)

s_x = []
s_y = []
t_xy = []

noe = len(edofs)
for i in range(noe):
    es, et = cfc.plants(ex[i,:], ey[i,:],ep,D,ed[i,:])
    s_x.append(es[0,0])
    s_y.append(es[0,1])
    t_xy.append(es[0,2])
# ------------------------------------------------------------------
# Vẽ sơ đồ biến dạng

etype = mesh.el_type
clmap="plasma"

cfv.figure()
cfv.draw_displacements(d*10,coords,edofs,ndof,etype,
                       draw_undisplaced_mesh=True,
                       title='Sơ đồ biến dạng')

# Vẽ biểu đồ ứng suất
cfv.figure()
cfv.draw_element_values(s_x,coords,edofs,ndof,etype, 
                        None, draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='s_x',
                        clmap=clmap)
cfv.axis('off')

cfv.figure()
cfv.draw_element_values(s_y,coords,edofs,ndof,etype, 
                        None,draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='s_y',
                        clmap=clmap)
cfv.axis('off')

cfv.figure()
cfv.draw_element_values(t_xy,coords,edofs,ndof,etype,
                        None,draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='t_xy',
                        clmap=clmap)
cfv.axis('off')
#cfv.colorbar()