from pprint import pprint
import csv
from math import sqrt, sin, cos, radians


cell_L_x = 62.840
cell_L_y = 62.840
cell_L_z = 62.840

cell_angle_a = radians(60)	
cell_angle_b = radians(60)
cell_angle_c = radians(60)

global ax, bx, by, cx, cy, cz
ax = cell_L_x
bx = cell_L_y*cos(cell_angle_c)
by = cell_L_y*sin(cell_angle_c)
cx = cell_L_z*cos(cell_angle_b)
cy = cell_L_z*(cos(cell_angle_a)-cos(cell_angle_b)*cos(cell_angle_c))/sin(cell_angle_c)
cz = sqrt(cell_L_z**2 - cx**2 - cy**2)



# mass cener of 4 points from the list 'pos'
def mass_center(pos,a,b,c,d):
    x = (pos[a][0] + pos[b][0] + pos[c][0] + pos[d][0]) / 4
    y = (pos[a][1] + pos[b][1] + pos[c][1] + pos[d][1]) / 4
    z = (pos[a][2] + pos[b][2] + pos[c][2] + pos[d][2]) / 4
    return x, y, z


# Judge the position is in or out of the cell
def InOut_Cell(pos):
    x = pos[0] 
    y = pos[1] 
    z = pos[2] 
    x_min = (bx*cz*y-(bx*cy-by*cx)*z)/by/cz
    x_max = x_min + ax
    y_min = z*cy/cz
    y_max = y_min + by
    z_min = 0
    z_max = cz

    if (x_min <= x < x_max) & (y_min <= y < y_max) & (z_min <= z < z_max):
        return True
    else:
        return False


def main():
    # read the position of O1 from cif file
    f = open('Framework_MIL101.pdb','r')
    posline = f.readlines()
    O1_list = []
    for line in posline:
        pos = line.split()
        if 'O1' in pos:
            xyz = [float(pos[4]), float(pos[5]), float(pos[6])]
            O1_list.append(xyz)
    f.close()
    num_O1 = len(O1_list)
    print('number of O1 in original cif file is :', num_O1)


    # expand periodic boundary
    for k in range(num_O1):
        for x in range(3):
            for y in range(3):
                for z in range(3):
                    if x==1 & y==1 & z==1:
                        continue
                    else:
                        # reference https://msestudent.com/triclinic-unit-cell/
                        posX = O1_list[k][0] + cx*(z-1) + bx*(y-1) + ax*(x-1)
                        posY = O1_list[k][1] + cy*(z-1) + by*(y-1)
                        posZ = O1_list[k][2] + cz*(z-1)
                        O1_list.append([posX, posY, posZ])

    num_O1 = len(O1_list)
    print('number of O1 in periodic cif file is :', num_O1)


    # write down the O1 position
    with open('O1_list.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(O1_list)

    # calculate the distance between each O1, close_O1 has the order number of O1 in the list of O1_list
    O1_distance = 15
    close_O1=[]
    for k in range(num_O1):
        l = []
        for i in range(num_O1):
            if k==i:
                continue
            r = sqrt((O1_list[k][0]-O1_list[i][0])**2 + (O1_list[k][1]-O1_list[i][1])**2 + (O1_list[k][2]-O1_list[i][2])**2)
            if r < O1_distance: # if two atoms is close enough
                l.append(i)
        close_O1.append(l)


    # find the tetrahedra strucuture
    tetra_list = []
    for p in range(num_O1): # select forst node
        for q in close_O1[p]: # select second node（you can make edge）
            com_node = list(set(close_O1[p]) & set(close_O1[q]))  # common node between 'p' node and 'q' node
            if len(com_node) == 2: # if there are two common node, you can make tetrahidra stracture
                tetra_list.append([p, q, com_node[0],com_node[1]])

    tetra_list = [sorted(l) for l in tetra_list]
    tetra_list = sorted(list(map(list, set(map(tuple, tetra_list))))) # erase the same structure
    print('number fo tetra structure in cif file is : ', len(tetra_list))


    # write down the tetra structure
    tetra_info = []
    numST  = 0
    for s in tetra_list:
        x, y, z = mass_center(O1_list, s[0], s[1], s[2], s[3])

        if InOut_Cell([x,y,z]): # if the point is inside of the cell
            numST += 1
            info = [x, y, z, numST, s[0], s[1], s[2], s[3]]
            tetra_info.append(info)


    with open('tetra_center.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(['center x', 'center y', 'center z', 'number', 'O1 num1','O1 num2', 'O1 num3', 'O1 num4'])
        writer.writerows(tetra_info)

    print('number of ST n the whole cell is :', len(tetra_info))



if __name__=='__main__':
    main()