import pandas as pd
import csv
from split import modify_PBDfile, read_line
from pprint import pprint

# judge whether it is in the sphere of the tetrahidra structure
def sphere_judge(r, xC,yC,zC, x_ST, y_ST, z_ST):
    dx = abs(xC - x_ST)
    dy = abs(yC - y_ST)
    dz = abs(zC - z_ST)

    # make the calcualtion Faster
    if (dx > r) or (dy > r) or (dz > r):
        return False

    dr = (dx**2 + dy**2 + dz**2)**(1/2)
    if dr <= r:
        return True
    else:
        return False


# count the histogram of the number of CO2 in the cell
def make_hist(counter):
    num_ST = len(counter)
    max_num_in_ST = 7
    hist = [0 for r in range(max_num_in_ST)]
    for c in counter:
        hist[c] += 1/num_ST*100 # normalized value in % representation
    return hist

    

def countST(P):
    # reading the position of the center of tetrahidra structure
    tetra = pd.DataFrame(pd.read_csv('tetra_center_PBC.csv'))

    # Radius of ST
    R = 4

    # read the division line of steps from pdb file
    model_line = read_line(fname = './movie_pdb/'+P+'.pdb')
    numSTEP = len(model_line) - 1
    print(' number of sampling gcmc steps : ',numSTEP)

    max_num_in_ST = 7
    hist_sum = [0 for r in range(max_num_in_ST)]

    for i in range(numSTEP):
        s = model_line[i] + 2
        row = model_line[i + 1] - model_line[i] - 3
        # print('read from line ', s+1, ' to line ', s+row)

        DataCo2 = modify_PBDfile(s, row, fname = './movie_pdb/'+P+'.pdb', step = i) # pass the data at step i
        count_CO2 = [0 for sts in range(0, 34)]

        # calcuate for each CO2 molecluaes in the cell 
        for index, posCO2 in DataCo2.iterrows():
            if (index-1)%3 == 0: # calculate only C atom in the trajectory
                x = float(posCO2[1])
                y = float(posCO2[2])
                z = float(posCO2[3])

                # iterate for each small cavity in the cell
                for k, posST in tetra.iterrows():
                    xST = float(posST['center x'])
                    yST = float(posST['center y'])
                    zST = float(posST['center z'])
                    numST = int(posST['number'])

                    # judge whether it is in the ST or not
                    if sphere_judge(R, x,y,z, xST, yST, zST):
                        count_CO2[int(numST)-1] += 1


        hist = (make_hist(count_CO2))
        print('Histogram of the number of CO2 at step ', str(i), ': ', hist)
        hist_sum = [hist_sum[r] + hist[r] for r in range(max_num_in_ST)]

    hist_sum = [hist_sum[r] / numSTEP for r in range(max_num_in_ST)]
    print('Histogram of the number of CO2 in the small cavity at pressure '+P+' Pa  : ', hist_sum)

    with open('count_CO2_in_ST.csv', 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow([P]+hist_sum)


if __name__ == '__main__':
    # target pressure
    Pressure = '1000000'
    countST(Pressure)