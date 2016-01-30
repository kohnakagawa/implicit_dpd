CXX = g++
#CXX = icpc

use_gpu_cuda = yes
#gpu_profile = yes

OPENMP = -DPARTICLE_SIMULATOR_THREAD_PARALLEL
ifeq ($(CXX),icpc)
OPENMP += -openmp 
endif
ifeq ($(CXX),g++)
OPENMP += -fopenmp
endif

WARNINGS = -Wall -Wunused-variable #-Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3 -DNDEBUG

INCLUDE = -I./FDPS/src

#COMMON_FLAGS = $(DEBUG)
COMMON_FLAGS = $(RELEASE)

ifeq ($(COMMON_FLAGS), $(DEBUG))
CUDA_DEBUG = -G
endif

COMMON_FLAGS += -std=c++11

ifeq ($(CXX),icpc)
OPT_FLAGS = -ipo -no-prec-div -xHOST
endif
ifeq ($(CXX),g++)
OPT_FLAGS = -ffast-math
endif

OBJECTS_DPD = main.o
OBJECTS_CMAKE = config_maker.o

ifeq ($(use_gpu_cuda),yes)
COMMON_FLAGS += -DENABLE_GPU_CUDA
CUDA_HOME = /usr/local/cuda
NVCC = $(CUDA_HOME)/bin/nvcc
NVCCFLAGS = $(COMMON_FLAGS) -ccbin=$(CXX) -arch=sm_35 -Xcompiler "$(COMMON_FLAGS) $(WARNINGS) $(OPT_FLAGS)" $(CUDA_DEBUG)

ifeq ($(gpu_profile),yes)
NVCCFLAGS += -lineinfo
endif

LIBRARY = -L$(CUDA_HOME)/lib64 -lcudart
INCLUDE += -I$(CUDA_HOME)/include/ -I$(CUDA_HOME)/samples/common/inc/
OBJECTS_DPD += f_calculator_gpu.o
endif

CXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(OPENMP) $(OPT_FLAGS)

TARGET_DPD = implicit_dpd.out
TARGET_CMAKE = config_maker.out

all:$(TARGET_DPD) $(TARGET_CMAKE)

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@
.SUFFIXES: .cu .o
.cu.o:
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@

$(TARGET_DPD): $(OBJECTS_DPD)
	$(CXX) $(CXX_FLAGS) $(OBJECTS_DPD) $(LIBRARY) -o $@

$(TARGET_CMAKE): $(OBJECTS_CMAKE)
	$(CXX) $(CXX_FLAGS) $(OBJECTS_CMAKE) $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS_DPD) $(OBJECTS_CMAKE) $(TARGET_DPD) $(TARGET_CMAKE) core.* *~
