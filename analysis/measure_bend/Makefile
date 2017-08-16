TARGET = bend.out

bend.out: measure_bend_coef.cpp
	g++ -Wall -Wextra -fopenmp -std=c++11 -O3 $< -o $@ -lfftw3

clean:
	rm -f $(TARGET) *~
