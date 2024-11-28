"""
Ví dụ 5.4
Tính vách đơn chịu tải trọng ngang
"""
# ------------------------------------------------------------------
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import calfem.geometry as cfg
import calfem.mesh as cfm
import numpy as np
# ------------------------------------------------------------------
# Tạo dạng hình học
g = cfg.Geometry()

bottom_edge = 20
left_edge = 21
top_point = 10

a,b = 1.2, 4.8
g.point([0, 0]) # point 0
g.point([a, 0]) # point 1
g.point([0, b]) # point 2
g.point([a, b], marker=top_point) # point 3

g.spline([0, 1], marker=bottom_edge)    # line 0
g.spline([2, 3])                        # line 1
g.spline([0, 2], marker=left_edge)      # line 2
g.spline([1, 3])                        # line 3

g.surface([0, 3, 1, 2])
g.c
cfv.draw_geometry(g)
# ------------------------------------------------------------------
# Tự động chia lưới
mesh = cfm.GmshMesh(g)
mesh.el_type = 2            # Element type
mesh.dofs_per_node = 2      # Degrees of freedom per node.
mesh.el_size_factor = 0.2   # Factor that changes element sizes.

coords, edofs, dofs, bdofs, elementmarkers = mesh.create()

ex, ey = cfc.coordxtr(edofs, coords, dofs)

# cfv.figure()
# cfv.draw_mesh(
#     coords=coords,
#     edof=edofs,
#     dofs_per_node=mesh.dofs_per_node,
#     el_type=mesh.el_type,
#     filled=True,
#     title='Lưới phần tử'
# )
# ------------------------------------------------------------------
# Số liệu tính toán
E,nu,t,q = 3e6, 0.25, 0.1, 10.
ptype = 1
ep = [ptype,t]
D = cfc.hooke(ptype,E,nu)

# lập ma traannj độ cứng của hệ
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
ndof = mesh.dofs_per_node
no_lpoint = len(bdofs[left_edge])/ndof
q_point = q*b/no_lpoint
cfu.applyforce(bdofs, f, left_edge, q_point, 1)

# Giải hệ phương trình
d, r = cfc.solveq(K,f,bc)

d_index = np.array(bdofs[top_point])-1
print('Chuyển vị ngang đỉnh vách',d[d_index][0])

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

cfv.figure()
cfv.draw_displacements(d*50,coords,edofs,ndof,etype,
                       draw_undisplaced_mesh=True,
                       title='Sơ đồ biến dạng')

# Vẽ biểu đồ ứng suất
cfv.figure()
cfv.draw_element_values(s_x,coords,edofs,ndof,etype, 
                        None, draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='s_x')
cfv.axis('off')

cfv.figure()
cfv.draw_element_values(s_y,coords,edofs,ndof,etype, 
                        None,draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='s_y')
cfv.axis('off')

cfv.figure()
cfv.draw_element_values(t_xy,coords,edofs,ndof,etype,
                        None,draw_elements=False, 
                        draw_undisplaced_mesh=False, 
                        title='t_xy')
cfv.axis('off')


