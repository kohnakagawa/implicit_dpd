CXX = g++
#CXX = icpc
NVCC = nvcc

ifeq ($(CXX),icpc)
OPENMP = -openmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif

WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual -Wunused-variable
DEBUG = -O0 -g -DDEBUG

ifeq ($(CXX),icpc)
RELEASE = -O3 -ipo -no-prec-div -xHOST
endif
ifeq ($(CXX),g++)
RELEASE = -O3 -march=native
endif

STDCPP11 = -std=c++11

#CFLAGS = $(DEBUG)
CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
#CFLAGS += $(OPENMP)
CFLAGS += $(STDCPP11)
CFLAGS += -DPAIRWISE_DPD

CUDA_HOME = /usr/local/cuda
INCLUDE = -isystem ../../FDPS/src -I../../include -I$(CUDA_HOME)/include/ -I$(CUDA_HOME)/samples/common/inc/

OBJECTS = prng_test.o

TARGET = prng_test.out saru_simd.out

all:$(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .o
prng_test.out: prng_test.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

saru_simd.out: saru_simd.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

cuda_prng.out: prng_test.cu
	$(NVCC) -O3 -std=c++11 -arch=sm_35 $(INCLUDE) $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~ cuda_prng.out
