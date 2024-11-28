"""
Ví dụ 5.2
Tạo lưới phần tử tam giác
"""

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

# ------------------------------------------------------------------
x1,y1 = 0,0
x2,y2 = 3,1
x3,y3 = 1.5,3
# ------------------------------------------------------------------
g = cfg.Geometry()
left_support = 10
right_support = 11
g.point([x1, y1], marker=left_support) 
g.point([x2, y2], marker=right_support) 
g.point([x3, y3])
g.spline([0, 1]) 
g.spline([1, 2], marker=20) 
g.spline([2, 0]) 
g.surface([0, 1, 2], marker=30)
cfv.figure()
cfv.draw_geometry(g)
# ------------------------------------------------------------------
mesh = cfm.GmshMesh(g)
mesh.el_type = 2            
mesh.dofs_per_node = 2      
mesh.el_size_factor = 0.5
coords,edofs,dofs,bdofs,elementmarkers = mesh.create()
cfv.figure()
cfv.draw_mesh(
    coords=coords,
    edof=edofs,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title='Lưới phần tử'
)

print(bdofs[left_support]) 
print(bdofs[right_support]) 

print(bdofs[20])