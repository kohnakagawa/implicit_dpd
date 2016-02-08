#!/bin/sh
export OMP_NUM_THREADS=4

mpiexec -np 4 ./mpi_bonded.out ./result_mpi/
./bonded.out ./result/

awk '{if(NR >= 3) print $0}' ./result_mpi/traject.xyz > tmp_mpi.txt
awk '{if(NR >= 3) print $0}' ./result/traject.xyz > tmp.txt
sort -n -k5 tmp_mpi.txt > sort.txt
if cmp -s file1 file2; then
   echo Success test
else
   echo Different
fi
rm tmp.txt tmp_mpi.txt sort.txt
