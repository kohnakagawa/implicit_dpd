#!/bin/sh

mpiexec -np 4 ./io_test_mpi.out ./result/

awk '{if(NR >= 3) print $0}' ./result/traject.xyz > tmp.txt
awk '{if(NR >= 3) print $0}' ./result/init_config.xyz > org.txt
sort -n -k5 tmp.txt > sort.txt
diff sort.txt org.txt
rm tmp.txt org.txt sort.txt

diff ./result/debug_dump_rank0.txt ./result/debug_dump_rank1.txt
diff ./result/debug_dump_rank1.txt ./result/debug_dump_rank2.txt
diff ./result/debug_dump_rank2.txt ./result/debug_dump_rank3.txt
