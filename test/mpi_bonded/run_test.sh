#!/bin/sh
export OMP_NUM_THREADS=4

#dir=random
#dir=flat
dir=cylind
#dir=sphere


mpiexec -np 4 ./mpi_bonded.out ./result_mpi/$dir
./bonded.out ./result_single/$dir

awk '{if(NR >= 3) print $0}' ./result_mpi/$dir/traject.xyz > tmp_mpi.txt
awk '{if(NR >= 3) print $0}' ./result_single/$dir/traject.xyz > tmp.txt
sort -n -k5 tmp_mpi.txt > sort.txt
if cmp -s sort.txt tmp.txt; then
   echo Success test.
else
   echo Different. Show maximum error.
fi

awk '{print $15, $16, $17}' sort.txt > 1.txt
awk '{print $15, $16, $17}' tmp.txt > 2.txt

R --vanilla --slave <<EOF
# check acc
print("acc error")
dat_mpi <- data.matrix(read.table("1.txt"))
dat     <- data.matrix(read.table("2.txt"))
diff    <- abs(dat_mpi - dat)
max(diff)

print("pressure error")
dat_mpi <- data.matrix(read.table("./result_mpi/$dir/pressure.txt"))
dat     <- data.matrix(read.table("./result/$dir/pressure.txt"))
diff    <- abs(dat_mpi - dat)
max(diff)
EOF

rm tmp_mpi.txt sort.txt tmp.txt 1.txt 2.txt