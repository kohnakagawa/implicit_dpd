#!/usr/bin/python

import os, sys, math

APL   = 0.10983521 # Area Per Lipid (chi 90 kappa 50)
Thick = 2.462681   # thickness

def main(dir_name):
    fin_a_name = os.path.join(dir_name, "elem_in_each.txt")
    fin_rg_name = os.path.join(dir_name, "rg.txt")
    fout_name = os.path.join(dir_name, "sc_ratio.txt")
    
    with open(fin_a_name, "r") as fin0, open(fin_rg_name, "r") as fin1, open(fout_name, "w") as fout:
        ade_data = fin0.readlines()[1:]
        rg_data  = fin1.readlines()
        num_line = len(rg_data)
        fout.write("# time ratio_to_rg ratio_to_actual_area")
        for i in range(num_line):
            time = float(rg_data[i].split()[0])
            rg   = float(rg_data[i].split()[1])
            num0 = int(ade_data[i].split()[1])
            num1 = int(ade_data[i].split()[2])
            
            delta_a0             = (num0 - num1) * APL
            area                 = (num0 + num1) * 0.5 * APL
            ratio_to_rg          = 2.0 * math.pi * Thick * rg / (delta_a0)
            ratio_to_actual_area = 2.0 * math.pi * Thick * math.sqrt(area / (4.0 * math.pi)) / delta_a0
            
            fout.write("%f %f %f\n" % (time, ratio_to_rg, ratio_to_actual_area))

if __name__ == "__main__":
    argvs = sys.argv
    argc  = len(argvs)
    if argc != 2:
        print "Usage"
        print "$ python %s target_dir\n" % argvs[0]
        quit()
    dir_name = argvs[1]
    main(dir_name)
