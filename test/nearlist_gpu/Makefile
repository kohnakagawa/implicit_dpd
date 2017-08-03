NVCC = nvcc

DEBUG = -O0 -g -DDEBUG

RELEASE = -O3 
STDCPP11 = -std=c++11

#CFLAGS = $(DEBUG)
CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(STDCPP11)
CFLAGS += -DPAIRWISE_DPD

CUDA_HOME = /usr/local/cuda

LIBRARY = -L$(CUDA_HOME)/lib64 -lcudart

INCLUDE = -isystem ../../FDPS/src -I../../include
INCLUDE += -I$(CUDA_HOME)/include/ -I$(CUDA_HOME)/samples/common/inc/ -Xcompiler "-O3"

TARGET = near_list.out

all:$(TARGET)

near_list.out: near_list.cu
	$(NVCC) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
