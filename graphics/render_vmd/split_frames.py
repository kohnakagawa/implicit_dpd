#!/usr/bin/python

import sys
import os.path

def writexyzform(fhand, lnum, dname, time, fid):
    fname = os.path.join(dname, "snapshot_id%d.xyz" % fid)
    with open(fname, "w") as fw:
        fw.write("%d\n" % lnum)
        fw.write("time %d\n" % time)
        for i in range(lnum):
            fw.write(fhand.readline())

def main(fin, dname):
    with open(fin, "r") as fhand:
        fid = 0
        while True:
            head = fhand.readline()
            if head == '':
                break
            lnum = int(head)
            time = int(fhand.readline().split()[1])
            writexyzform(fhand, lnum, dname, time, fid)
            fid += 1
        print fid

if __name__ == "__main__":
    argvs = sys.argv
    argc  = len(argvs)
    
    if argc != 2:
        print "Please specifiy a target directory name"
        print "argv[1] is target directory name"
        quit()
    
    fin = os.path.join(argvs[1], "traject.xyz")
    if not os.path.exists(fin):
        print "Cannot find tranjectory file in %s" % fin
        quit()
    
    main(fin, argvs[1])
