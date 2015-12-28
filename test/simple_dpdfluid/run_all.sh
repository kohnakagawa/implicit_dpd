#!/bin/sh

rm -f rho_vs_press.txt kin_temp.txt

for i in "1.0" "2.0" "3.0" "4.0" "5.0" "6.0" "7.0"
do
    ./groot_warren.out ./test_result/rho$i 200000 pressure &
done

for i in "0.01" "0.005" "0.001"
do
    ./groot_warren.out ./test_result/dt$i 200000 temperature &
done

wait

export OMP_NUM_THREADS=4

for i in "1.0" "2.0" "3.0" "4.0" "5.0"
do
    ./groot_warren_omp.out ./test_result/rho$i 200000 pressure &
done

wait

for i in "6.0" "7.0"
do
    ./groot_warren_omp.out ./test_result/rho$i 200000 pressure &
done

for i in "0.01" "0.005" "0.001"
do
    ./groot_warren_omp.out ./test_result/dt$i 200000 temperature &
done

wait