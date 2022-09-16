import pandas as pd
from math import sqrt
from FindTetraInCell import ax, bx, by, cx, cy, cz


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
                    pos = pos.append(imT)
    return pos


or_Tetra = pd.DataFrame(pd.read_csv('tetra_center.csv'))
ex_Tetra = repro_pbd(or_Tetra)
ex_Tetra.to_csv('tetra_center_PBC.csv')