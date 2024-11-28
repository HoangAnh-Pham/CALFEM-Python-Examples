"""
Ví dụ 5.3
Tính vách tam giác và khảo sát hội tụ của kết quả
"""
# ------------------------------------------------------------------
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import calfem.geometry as cfg
import calfem.mesh as cfm
import numpy as np
# ------------------------------------------------------------------
x1,y1 = 0,0
x2,y2 = 3,1
x3,y3 = 1.5,3
# ------------------------------------------------------------------
# Tạo dạng hình học
g = cfg.Geometry()
left_support = 10
right_support = 11
top_point = 12

g.point([x1, y1], marker=left_support) 
g.point([x2, y2], marker=right_support) 
g.point([x3, y3], marker=top_point)

g.spline([0, 1]) 
g.spline([1, 2]) 
g.spline([2, 0]) 

g.surface([0, 1, 2])

# cfv.draw_geometry(g)
# ------------------------------------------------------------------
# Tự động chia lưới
mesh = cfm.GmshMesh(g)
mesh.el_type = 2            
mesh.dofs_per_node = 2      
mesh.el_size_factor = 0.02
coords,edofs,dofs,bdofs,elementmarkers = mesh.create()

ex, ey = cfc.coordxtr(edofs,coords,dofs)

cfv.figure()
cfv.draw_mesh(
    coords=coords,
    edof=edofs,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Lưới phần tử "
)
# ------------------------------------------------------------------
# Tính toán
E,nu,t,bx,by = 3e6, 0.3, 0.02, 5e3, 5e3
ptype = 1
ep = [ptype,t]
eq = [[bx],[by]]
D = cfc.hooke(ptype,E,nu)

nod = np.size(dofs)
K = np.zeros([nod, nod])
f = np.zeros([nod,1])
for eldof, elx, ely in zip(edofs, ex, ey):
    Ke,fe = cfc.plante(elx, ely, ep, D, eq)
    cfc.assem(eldof, K, Ke, f, fe)
# ------------------------------------------------------------------
bc = np.array([],'i')
bcVal = np.array([],'f')    
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 
                        left_support, 0.0, 0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 
                        right_support, 0.0, 2)
print(bc)
print(bcVal)
# ------------------------------------------------------------------
d, r = cfc.solveq(K,f,bc)

d_index = np.array(bdofs[top_point])-1
print(d[d_index])
print(len(edofs))




