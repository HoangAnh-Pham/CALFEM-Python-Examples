'''
Ví dụ 6.1
Chương trình tạo lưới phần tử chữ nhật
'''

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

def rectangular_mesh(a, b, m, n):

    # Define the number of nodes in x and y directions
    num_nodes_x = m + 1
    num_nodes_y = n + 1

    # Define the coordinates of the nodes
    x_coords = np.linspace(0, a, num_nodes_x)
    y_coords = np.linspace(0, b, num_nodes_y)
    coords = np.array([(x, y) 
                       for y in y_coords
                       for x in x_coords])

    # ----------------------------------------------------------------
    ndof = 3
    nnode = num_nodes_x*num_nodes_y
    dofs = np.arange(nnode*ndof) + 1
    dofs = np.array(dofs.reshape(nnode,ndof))
    # ----------------------------------------------------------------
    # Define the connectivity of the elements
    edofs = []
    for j in range(n):
        for i in range(m):
            n1 = j * num_nodes_x + i
            n2 = n1 + 1
            n3 = n1 + num_nodes_x
            n4 = n3 + 1
            dof = np.array([dofs[n1],
                            dofs[n2],
                            dofs[n4],
                            dofs[n3]])
            edofs.append(dof.flatten())
    edofs = np.array(edofs)
    # ----------------------------------------------------------------
    # Plot the mesh
    ex, ey = cfc.coordxtr(edofs, coords, dofs)
    cfv.eldraw2(ex, ey)
    # ----------------------------------------------------------------
    # Boundary dofs
    b1 = [i 
          for i in range(num_nodes_x)] # bottom
    b2 = [num_nodes_x*n + i 
          for i in range(num_nodes_x)] # top
    b3 = [j*num_nodes_x 
          for j in range(num_nodes_y)] # left
    b4 = [(j+1)*num_nodes_x - 1 
          for j in range(num_nodes_y)] # right

    bdof1 = []
    for node in b1:
        bdof1.append(dofs[node])
    bdof1 = np.array(bdof1)
    bdof1 = bdof1.flatten()

    bdof2 = []
    for node in b2:
        bdof2.append(dofs[node])
    bdof2 = np.array(bdof2)
    bdof2 = bdof2.flatten()

    bdof3 = []
    for node in b3:
        bdof3.append(dofs[node])
    bdof3 = np.array(bdof3)
    bdof3 = bdof3.flatten()

    bdof4 = []
    for node in b4:
        bdof4.append(dofs[node])
    bdof4 = np.array(bdof4)
    bdof4 = bdof4.flatten()

    bdofs = {}
    bdofs[1] = bdof1
    bdofs[2] = bdof2
    bdofs[3] = bdof3
    bdofs[4] = bdof4
    
    return coords, edofs, dofs, bdofs
