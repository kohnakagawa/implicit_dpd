#!/usr/bin/env python

import glob
import os
import sys
import re
import shutil

def copy_with_log(src, dst):
    print "%s -> %s" % (src, dst)
    shutil.copy(src, dst)

def get_backup_number(root_dir, f_back_pattern):
    backup_nums = set()
    for f in os.listdir(root_dir):
        f_number = re.search(f_back_pattern, f)
        if f_number is None:
            continue
        matched_num = f_number.group(3)
        backup_nums.add(matched_num)

    if len(backup_nums) == 0:
        return 0
    else:
        return int(max(backup_nums)) + 1

def make_backup(root_dir):
    f_back_pattern = r'(\w+)\.(\w+)\.(\d+)'

    num = get_backup_number(root_dir, f_back_pattern)
    for f in os.listdir(root_dir):
        sim_data = re.search(f_back_pattern, f)
        if sim_data is not None: # skip backup files
            continue
        f_back = f + "." + str(num)

        f = os.path.join(root_dir, f)
        f_back = os.path.join(root_dir, f_back)
        copy_with_log(f, f_back)

def make_init_config(root_dir):
    init_config = os.path.join(root_dir, "init_config.xyz")
    fin_config  = os.path.join(root_dir, "fin_config.xyz")
    if os.path.getsize(fin_config) == 0:
        print "WARNING! there is no trajectory data in %s" % fin_config
        sys.exit(1)
    copy_with_log(fin_config, init_config)

def main(argv):
    if len(argv) != 2:
        print "Usage: %s root_dir" % argv[0]
        sys.exit(1)
    root_dir = argv[1]
    make_backup(root_dir)
    make_init_config(root_dir)

if __name__ == "__main__":
    main(sys.argv)
