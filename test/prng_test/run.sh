#!/bin/sh

make
./prng_test.out
R --vanilla < show_dist.R
evince Rplots.pdf