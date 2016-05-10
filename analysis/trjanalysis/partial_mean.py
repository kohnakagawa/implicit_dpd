#!/usr/bin/python

import sys
import os
import numpy as np

def main(dname, tar_id, min_z, max_z):
    fin  = os.path.join(dname, "slice_rad.txt")
    fout_name = "mean_tube_rad_block%d.txt" % tar_id
    fout = os.path.join(dname, fout_name)
    with open(fin, "r") as f_red, open(fout, "w") as f_wrt:
        block_id = 0
        trad  = np.array([])
        zcord = np.array([])
        for line in f_red.readlines():
            if (block_id == tar_id) & (line != '\n'):
                line_elem = line.split()
                zcord = np.append(zcord, float(line_elem[0]))
                trad  = np.append(trad,  float(line_elem[1]))
            if line == '\n':
                block_id += 1
        tube_rad = np.mean(trad[(zcord > min_z) & (zcord < max_z)])
        f_wrt.write("Tube radius is %.10g\n" % tube_rad)
        print "Tube radius is %.10g" % tube_rad

if __name__  == "__main__":
    argvs = sys.argv
    argc  = len(argvs)
    if argc != 5:
        print "Usage: $ python %s target_dir block_id min_z max_z" % argvs[0]
        quit()
    main(argvs[1], int(argvs[2]), float(argvs[3]), float(argvs[4]))
