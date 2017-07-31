#!/usr/bin/env python
import os
import sys

def stringize_value(val):
    if isinstance(val, list):
        return " ".join(map(lambda x: str(x), val))
    else:
        return str(val)

def gen_input_script(param, root_dir, name):
    fname = os.path.join(root_dir, name)
    with open(fname, "w") as f:
        for key, value in param.items():
            f.write(" ".join([key, stringize_value(value)]) + "\n")

def run_input_base():
    return {
        "box_leng"     : [0.0, 0.0, 0.0],
        "init_amp_num" : 0,
        "sol_num"      : 0,
        "dt"           : 0.005,
        "cf_r"         : [3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
        "cf_s"         : 19.0,
        "cf_b"         : 5.0,
        "chi"          : 30.0,
        "kappa"        : 100.0,
        "rho_co"       : 40.0,
        "all_time"     : 500000,
        "step_mic"     : 50000,
        "step_mac"     : 10,
        "p_thresld"    : 0.0,
        "eps"          : 0.0,
        "max_amp_num"  : 0,
        "influ_hei"    : 0.0,
        "influ_dep"    : 0.0,
        "influ_vol"    : 0.0,
        "upside_rad"   : 0.0,
        "influ_rad"    : 0.0,
        "influ_grd"    : 0.0,
        "core_amp_id"  : [0, 1, 2],
    }

def conf_input_base():
    return {
        "Temperature"  : 1.0,
        "box_leng"     : [0.0, 0.0, 0.0],
        "amp_num"      : 0,
        "mode"         : "flat",
        "head_unit"    : 4,
        "tail_unit"    : 12,
        "bond_leng"    : 0.5,
        "sph_rad"      : 0.0,
        "upper_the"    : 0.0,
        "cyl_l"        : 0.0,
        "cyl_r"        : 0.0,
        "in_out_rat"   : 0.0,
        "sol_num"      : 0,
        "lip_len"      : 0.5,
    }

def generate_inputs(root_dir):
    with_solvent = False
    leng = [27.0, 27.5, 28.0, 28.5, 29.0]
    area_per_lipid = 0.333

    dir_names = ["box" + str(l) + "x" + str(l) for l in leng]

    for l, d_name in zip(leng, dir_names):
        run_params  = run_input_base()
        conf_params = conf_input_base()
        if with_solvent:
            conf_params["with_solvent"] = True

        box_leng  = [l, l, l]
        amp_num   = int(2 * (l * l / area_per_lipid))

        run_params["box_leng"] = box_leng
        run_params["init_amp_num"] = amp_num

        conf_params["box_leng"] = box_leng
        conf_params["amp_num"]  = amp_num

        param_dir = os.path.join(root_dir, d_name)
        if not os.path.exists(param_dir):
            os.mkdir(param_dir)

        gen_input_script(run_params, param_dir, "run_param.txt")
        gen_input_script(conf_params, param_dir, "config_param.txt")

def main(argv):
    if len(argv) != 2:
        print "Usage: %s root_dir" % argv[0]
        sys.exit(1)
    root_dir = argv[1]
    generate_inputs(root_dir)

if __name__ == "__main__":
    main(sys.argv)
