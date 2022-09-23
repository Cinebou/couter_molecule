import pandas as pd
from math import sqrt
from FindTetraInCell import ax, bx, by, cx, cy, cz


global surface_depth
# depth of the surrounding surface of the original cell
surface_depth = 5 

def judge_surface(x, y, z):
    # surface depth is extended in the virtical direction of the wall 
    x_min_surface = (bx*cz*y-(bx*cy-by*cx)*z)/by/cz - surface_depth*sqrt((by*cz)**2+(bx*cz)**2+(bx*cy-by*cx)**2)/by/cz
    x_max_surface = x_min_surface + ax + surface_depth*sqrt((by*cz)**2+(bx*cz)**2+(bx*cy-by*cx)**2)/by/cz + 3
     
    y_min_surface = z*cy/cz - surface_depth*sqrt(cy**2+cz**2)/cy +3
    y_max_surface = y_min_surface + by + surface_depth*sqrt(cy**2+cz**2)/cy

    # z surface is pararell to the original xyz coordinates 
    z_min_surface = 0 - surface_depth
    z_max_surface = cz + surface_depth

    if (x_min_surface <= x <= x_max_surface) & (y_min_surface <= y <= y_max_surface) & (z_min_surface <= z <= z_max_surface):
        return True
    else:
        return False



def repro_pbd(pos):
    for k in range(len(pos)):
        for x in range(3):
            for y in range(3):
                for z in range(3):
                    if x==1 & y==1 & z==1:
                        continue
                    imT = pos.iloc[k].copy(deep=True)
                    imT['center x'] += cx*(z-1) + bx*(y-1) + ax*(x-1)
                    imT['center y'] += cy*(z-1) + by*(y-1)
                    imT['center z'] += cz*(z-1)

                    # only the surface of the surrounding cell can be listed
                    if judge_surface(imT['center x'], imT['center y'], imT['center z']):
                        pos = pos.append(imT)
    return pos


def main():
    or_Tetra = pd.DataFrame(pd.read_csv('tetra_center.csv'))
    ex_Tetra = repro_pbd(or_Tetra)
    print('number of small cavity in the file is  : ', len(ex_Tetra))
    ex_Tetra.to_csv('tetra_center_PBC.csv')


if __name__=='__main__':
    main()