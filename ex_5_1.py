
"""
Ví dụ 5.1
Tính tấm phẳng tam giác bằng một phần tử
"""
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np

E,nu,t,bx,by = 3e6, 0.3, 0.02, 5e3, 5e3

edof = np.array([1,2,3,4,5,6])

x1,y1 = 0,0
x2,y2 = 3,1
x3,y3 = 1.5,3

ex = [x1,x2,x3]
ey = [y1,y2,y3]
ep = [1,t]
eq = [[bx],[by]]

D = cfc.hooke(1,E,nu)

Ke,fe = cfc.plante(ex,ey,ep,D,eq)
cfu.disp_h1('Ma trận độ cứng')
cfu.disp_array(Ke)
cfu.disp_h1('Véc tơ lực nút')
cfu.disp_array(fe)

K = np.matrix(np.zeros((6,6)))
f = np.matrix(np.zeros((6,1)))
cfc.assem(edof,K,Ke,f,fe)

bc = np.array([1,2,4])

d,r = cfc.solveq(K,f,bc)

ed = cfc.extract_eldisp(edof,d)

es,et = cfc.plants(ex,ey,ep,D,ed)
ef = cfc.plantf(ex,ey,ep,es)

ndis = d.reshape(3,2).copy()
cfu.disp_array(ndis,['u_x','u_y'])
cfu.disp_array(es,['s_x','s_y','t_xy'])
cfu.disp_array(np.array([ef]),
               ['f1','f2','f3','f4','f5','f6'])


