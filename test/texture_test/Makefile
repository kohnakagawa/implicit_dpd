CXX = g++
CUDA_HOME = /usr/local/cuda
NVCC = $(CUDA_HOME)/bin/nvcc
NVCCFLAGS = -O0 -g -G -ccbin=$(CXX) -arch=sm_35 -Xcompiler "-g -O0" -std=c++11

all:test.out

test.out: texture_test.cu
	$(NVCC) $(NVCCFLAGS) $< -o $@

clean:
	rm -f test.out core.* *~
