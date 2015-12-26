#!/bin/sh

for i in "1.0" "2.0" "3.0" "4.0" "5.0" "6.0" "7.0"
do
    ./groot_warren.out ./result/rho$i 200000 pressure &
done

for i in "0.01" "0.005" "0.001"
do
    ./groot_warren.out ./result/dt$i 200000 temperature &
done

wait