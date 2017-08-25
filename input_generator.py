#!/usr/bin/env python
import os
import sys
import math
import subprocess

def stringize_value(val):
    if isinstance(val, list):
        return " ".join(map(lambda x: str(x), val))
    else:
        return str(val)

def write_script(param, root_dir, name):
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
        "cf_b_rigid"   : 5.0,
        "vbb"          : 0.1,
        "chi"          : 30.0,
        "kappa"        : 100.0,
        "rho_co"       : 40.0,
        "all_time"     : 500000,
        "step_mic"     : 50000,
        "step_mac"     : 1,
        "p_thresld"    : -0.2,
        "eps"          : 0.2,
        "max_amp_num"  : 0,
        "influ_hei"    : 0.8,
        "influ_dep"    : 4.0,
        "influ_vol"    : 162.0,
        "upside_rad"   : 3.9,
        "influ_rad"    : 2.0,
        "influ_grd"    : -0.5,
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

def gen_inputs(root_dir):
    with_solvent = False
    leng = [27.0, 27.5, 28.0, 28.5, 29.0]
    area_per_lipid = 0.333

    dir_names = [os.path.join(root_dir, "box" + str(l) + "x" + str(l)) for l in leng]
    for l, param_dir in zip(leng, dir_names):
        run_params  = run_input_base()
        conf_params = conf_input_base()
        if with_solvent:
            conf_params["with_solvent"] = True

        box_leng  = [l, l, l]
        amp_num   = int(2 * (l * l / area_per_lipid))

        run_params["box_leng"] = box_leng
        run_params["init_amp_num"] = amp_num
        run_params["max_amp_num"] = amp_num

        conf_params["box_leng"] = box_leng
        conf_params["amp_num"]  = amp_num

        if not os.path.exists(param_dir):
            os.mkdir(param_dir)

        write_script(run_params, param_dir, "run_param.txt")
        write_script(conf_params, param_dir, "config_param.txt")
    return dir_names

def get_submit_command(scheduler):
    if scheduler == "slurm":
        return "sbatch"
    elif scheduler == "pbs":
        return "xqsub"
    else:
        return "nohup"

def submit_template(scheduler, script, prog_name, prog_arg):
    cmd = get_submit_command(scheduler)
    prog_arg = stringize_value(prog_arg)
    if scheduler == "slurm":
        return "%s %s %s %s\n" % (cmd, script, prog_name, prog_arg)
    elif scheduler == "pbs":
        return "%s -sargs \"%s %s\" %s\n" % (cmd, prog_name, prog_arg, script)
    else:
        return "%s %s %s %s &\n" % (cmd, script, prog_name, prog_arg)

def gen_submit_script(prog_name, prog_args, scheduler, run_script, submit_script):
    if prog_name[0:2] != './':
        prog_name = "./" + prog_name

    with open(submit_script, "w") as f:
        f.write("#!/bin/sh\n")
        if scheduler != "pbs_bulk":
            for prog_arg in prog_args:
                f.write(submit_template(scheduler, run_script, prog_name, prog_arg))
        else:
            f.write("qsub %s\n" % run_script)

def gen_run_script(scheduler, nprocs, nthreads, nnodes, qname, run_script, joblist_file):
    with open(run_script, "w") as f:
        f.write("#!/bin/sh\n")

        if scheduler == "slurm":
            check_qname(qname, nox_qnames())
            options = slurm_options(nprocs, nthreads, nnodes, qname)
            write_slurm_options(options, f)
        elif scheduler == "pbs":
            check_qname(qname, sekirei_qnames())
            options = qsub_options(nprocs, nthreads, nnodes, qname)
            write_qsub_options(options, f)
            walltime = get_walltime_from_queue(qname)
            jobname = "Test"
            options = pbs_options(walltime, jobname)
            write_pbs_options(options, f)
        elif scheduler == "pbs_bulk":
            check_qname(qname, sekirei_qnames())
            options = qsub_bulk_options(nnodes, qname)
            write_qsub_options(options, f)
            walltime = get_walltime_from_queue(qname)
            jobname = "Bulk"
            options = pbs_options(walltime, jobname)
            write_pbs_options(options, f)
        else:
            print "job scheduler is not specified"

        if nthreads != 1:
            f.write("export OMP_NUM_THREADS=%d\n" % nthreads)

        if nprocs != 1 and scheduler != "pbs_bulk":
            if scheduler == "pbs":
                f.write("mpijob ")
            else:
                f.write("mpirun -np %d " % nprocs)

        if scheduler == "pbs":
            f.write("${SCRIPT_ARGS}\n")
        elif scheduler == "pbs_bulk":
            f.write("bulkjob %s\n" % joblist_file)
        else:
            f.write("$@\n")

def get_program_type(nprocs, nthreads):
    if nprocs != 1 and nthreads != 1:
        return 'H'
    elif nprocs == 1:
        return 'O'
    elif nthreads == 1:
        return 'M'
    else:
        return 'S'

def get_parallel_procs(prog_type, nprocs, nthreads):
    if prog_type == 'H':
        return 'H' + str(nthreads)
    elif prog_type == 'O':
        return nthreads
    elif prog_type == 'M':
        return nprocs
    elif prog_type == 'S':
        return 1
    else:
        print "Unknown program type %s\n" % prog_type
        sys.exit()

def bulkjob_list_template(prog_name, prog_arg, nprocs, nthreads, exec_dir):
    prog_name = os.path.abspath(prog_name)
    exec_dir  = os.path.abspath(exec_dir)
    command_line = "%s %s" % (prog_name, prog_arg)
    program_type = get_program_type(nprocs, nthreads)
    parallel_procs = str(get_parallel_procs(program_type, nprocs, nthreads))
    return ";".join((command_line, parallel_procs, program_type, exec_dir)) + "\n"

def gen_bulkjob_list(prog_name, prog_args, nprocs, nthreads, joblist):
    with open(joblist, "w") as f:
        for prog_arg in prog_args:
            f.write(bulkjob_list_template(prog_name, \
                                          stringize_value(prog_arg), \
                                          nprocs, nthreads, "./"))

def check_qname(name, qname_sets):
    if not name in qname_sets:
        print "%s unknown name" % name
        sys.exit()

def nox_qnames():
    return ("defq", "hwq", "sbq", "wmq")

def sekirei_qnames():
    return sekirei_qnames_F() + sekirei_qnames_L()

def sekirei_qnames_F():
    return ("F4cpu", "F36cpu", "F144cpu", "F18acc", "F72acc", "F2fat")

def sekirei_qnames_L():
    return ("L4cpu", "L36cpu", "L144cpu", "L18acc", "L2fat")

def qsub_options(nprocs, nthreads, nnodes, queue):
    return {
        "-queue" : queue,
        "-node"  : nnodes,
        "-mpi"   : nprocs,
        "-omp"   : nthreads,
        "-place" : "distribute",
        "-over"  : "false",
    }

def qsub_bulk_options(nnodes, queue):
    return {
        "-queue" : queue,
        "-node"  : nnodes,
    }

def pbs_options(walltime, jobname):
    return {
        "-l" : walltime,
        "-N" : jobname
    }

def slurm_options(nprocs, nthreads, nnodes, queue):
    return {
        "-p" : queue,
        "-n" : nprocs,
        "-N" : nnodes,
        "-c" : nthreads,
        "-o" : "./out_%J.log",
        "-e" : "./err_%J.log",
    }

def get_walltime_from_queue(qname):
    walltime="walltime="
    if qname in sekirei_qnames_F():
        walltime += "24:00:00"
    elif qname in sekirei_qnames_L():
        walltime += "120:00:00"
    else:
        print "%s unknown queue name"
        sys.exit(1)
    return walltime

def write_jobscheduler_options(options, f_hand, prefix):
    for key, value in options.items():
        f_hand.write(" ".join([prefix, key, stringize_value(value)]) + "\n")

def write_slurm_options(options, f_hand):
    return write_jobscheduler_options(options, f_hand, "#SBATCH")

def write_qsub_options(options, f_hand):
    return write_jobscheduler_options(options, f_hand, "#QSUB")

def write_pbs_options(options, f_hand):
    return write_jobscheduler_options(options, f_hand, "#PBS")

def omp_mpi_node_info():
    return {
        "nprocs"   : 1,
        "nthreads" : 4,
        "nnodes"   : 1,
    }

def gen_scripts(root_dir, dir_names, prog_name):
    run_script    = os.path.join(root_dir, "run.sh")
    submit_script = os.path.join(root_dir, "submit.sh")
    joblist_file  = os.path.join(root_dir, "joblist.txt")

    # for nox00
    # scheduler = "slurm"
    # qname     = "hwq"

    # for sekirei no bulk job
    # scheduler = "pbs"
    # qname     = "F18acc"

    # for sekirei bulk job
    scheduler = "pbs_bulk"
    qname     = "F36cpu"

    # other
    # scheduler = None
    # qname     = None

    # program arguments
    observe_beg_time = 100000
    prog_args = [[d, observe_beg_time] for d in dir_names]

    # omp mpi node info
    run_info = omp_mpi_node_info()

    gen_submit_script(prog_name = prog_name,        \
                      prog_args = prog_args,        \
                      scheduler = scheduler,        \
                      run_script = run_script,      \
                      submit_script = submit_script)
    gen_run_script(scheduler = scheduler,            \
                   nprocs    = run_info["nprocs"],   \
                   nthreads  = run_info["nthreads"], \
                   nnodes    = run_info["nnodes"],   \
                   qname     = qname,                \
                   run_script = run_script,          \
                   joblist_file = joblist_file)
    if scheduler == "pbs_bulk":
        gen_bulkjob_list(prog_name = prog_name,           \
                         prog_args = prog_args,           \
                         nprocs = run_info["nprocs"],     \
                         nthreads = run_info["nthreads"], \
                         joblist = joblist_file)

    os.chmod(run_script, 0775)
    os.chmod(submit_script, 0775)
    if scheduler == "pbs_bulk":
        os.chmod(joblist_file, 0775)

def gen_initconfig(confmaker_path, dir_names):
    confmaker = os.path.join(confmaker_path, "config_maker.out")
    if not os.path.exists(confmaker):
        print "%s does not exist" % confmaker
        print "init_config.xyz is not generated."
        print "You should specify confmaker_path manually."
        return

    is_first   = True
    be_written = True
    for d in dir_names:
        config = os.path.join(d, "init_config.xyz")
        if os.path.exists(config) and is_first:
            is_first   = False
            be_written = False
            print "init_config.xyz exits at %s" % d
            print "Do you overwrite this file [y/N]?"
            y_or_n = raw_input()
            if y_or_n == "y":
                be_written = True

        if be_written:
            cmd = [confmaker, d]
            subprocess.call(cmd)

def main(argv):
    if len(argv) != 4:
        print "Usage: %s root_dir confmaker_dir program_name" % argv[0]
        sys.exit(1)
    root_dir        = argv[1]
    confmaker_dir   = argv[2]
    program_name    = argv[3]
    input_dir_names = gen_inputs(root_dir)
    gen_scripts(root_dir, input_dir_names, program_name)
    gen_initconfig(confmaker_dir, input_dir_names)

if __name__ == "__main__":
    main(sys.argv)
