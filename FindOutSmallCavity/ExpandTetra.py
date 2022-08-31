import pandas as pd
from math import sqrt

global cell_L
cell_L = 62.8339

def repro_pbd(pos):
    for k in range(len(pos)):
        for x in range(3):
            for y in range(3):
                for z in range(3):
                    if x==1 & y==1 & z==1:
                        continue
                    imT = pos.iloc[k].copy(deep=True)
                    imT['center x'] += cell_L * ((x-1) + (y-1)/2 + (z-1)/2)
                    imT['center y'] += cell_L * ((y-1)*sqrt(3)/2 + (z-1)/2/sqrt(3))
                    imT['center z'] += cell_L * (z-1)*sqrt(2/3)
                    pos = pos.append(imT)
                    

    return pos


or_Tetra = pd.DataFrame(pd.read_csv('tetra_center.csv'))
ex_Tetra = repro_pbd(or_Tetra)
print(ex_Tetra)
ex_Tetra.to_csv('tetra_center_PBC.csv')