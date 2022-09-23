import pandas as pd
from pprint import pprint

def main():
    cav_pos = pd.read_csv('tetra_center_PBC.csv')
    tcl_lines = ''
    for index, row in cav_pos.iterrows():
        tcl_lines += 'draw sphere {{{} {} {}}} radius 4\n'.format(row['center x'], row['center y'], row['center z'])

    with open('small_cavity_sphere.tcl', 'w') as f:
        f.write(tcl_lines)
    return 0


if __name__=='__main__':
    main()