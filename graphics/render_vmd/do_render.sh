#!/bin/sh

if [ $# -ne 1 ]; then
  echo "please specify target directory"
  exit 1
fi

time=$(python ./split_frames.py $1)

vmd -e render_eps.tcl

clean up
rm -f snapshot_id*.xyz