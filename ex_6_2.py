'''
Ví dụ 6.2
Tính độ võng của tấm chữ nhật 
'''

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv
import rectangular_mesh as rm

# Only es is included as it's the only needed output
def platrs(ex,ey,ep,D,ed):
    
    Lx=ex[2]-ex[0]
    Ly=ey[2]-ey[0]
    t=ep[0]

    D=((t**3)/12)*D

    A1=D[1,1]/2/Ly
    A2=D[0,0]/2/Lx
    A3=D[0,1]/2/Ly
    A4=D[0,1]/2/Lx
    A5=D[2,2]/2/Ly
    A6=D[2,2]/2/Lx
    A7=4*D[2,2]/Lx/Ly

    B1=6*D[1,1]/Ly/Ly/Ly
    B2=6*D[0,0]/Lx/Lx/Lx
    B3=-3*D[1,1]/Ly/Ly
    B4=3*D[0,0]/Lx/Lx
    B5=D[0,1]/Lx/Ly

    mx=A3*(-ed[1]-ed[4]+ed[7]+ed[10])+A2*(ed[2]-ed[5]-ed[8]+ed[11])
    my=A1*(-ed[1]-ed[4]+ed[7]+ed[10])+A4*(ed[2]-ed[5]-ed[8]+ed[11])
    mxy=A6*(ed[1]-ed[4]-ed[7]+ed[10])+A5*(-ed[2]-ed[5]+ed[8]+ed[11])+A7*(ed[0]-ed[3]+ed[6]-ed[9])

    m1=0.5*(mx+my)+np.sqrt(0.25*(mx-my)**2+mxy**2)
    m2=0.5*(mx+my)-np.sqrt(0.25*(mx-my)**2+mxy**2)
    alfa=0.5*180/np.pi*np.arctan2(mxy,(mx-my)/2)

    vx=B5*(-ed[1]+ed[4]-ed[7]+ed[10])+B4*(ed[2]+ed[5]+ed[8]+ed[11])+B2*(-ed[0]+ed[3]+ed[6]-ed[9])
    vy=B3*(ed[1]+ed[4]+ed[7]+ed[10])+B5*(ed[2]-ed[5]+ed[8]-ed[11])+B1*(-ed[0]-ed[3]+ed[6]+ed[9])

    es=np.transpose(np.array([mx, my, mxy, vx, vy]))

    return es

# -----------------------------------------------------------------
# Ví dụ 6.2
# -----------------------------------------------------------------
a, b, m, n = 4, 3, 24, 24
coords,edofs,dofs,bdofs = rm.rectangular_mesh(a,b,m,n)
ex, ey = cfc.coordxtr(edofs, coords, dofs)
# -----------------------------------------------------------------
E,nu,t,q = 3e6, 0.3, 0.1, 10
ptype = 1
ep = [t]
eq = q
D = cfc.hooke(ptype,E,nu)
# -----------------------------------------------------------------
nod = np.size(dofs)
K = np.zeros([nod,nod])
f = np.zeros([nod,1])

Ke, fe = cfc.platre(ex[0], ey[0], ep, D, eq)
fe = fe.T
cfc.assem(edofs, K, Ke, f, fe)
# -----------------------------------------------------------------    
bc1 = bdofs[1][0::3]
bc2 = bdofs[2][0::3]
bc3 = bdofs[3][0::3]
bc4 = bdofs[4][0::3]
bc = np.array([bc1,bc2,bc3,bc4])
bc = bc.flatten()

d, r = cfc.solveq(K,f,bc)
w = d.reshape(len(dofs),3)[:,0]
w_max = max(abs(w))
print(w_max)
