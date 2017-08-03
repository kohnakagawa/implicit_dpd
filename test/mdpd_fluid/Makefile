CXX = g++
# CXX = icpc

ifeq ($(CXX),icpc)
OPENMP = -openmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif

WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual -Wunused-variable
DEBUG = -O0 -g -DDEBUG

ifeq ($(CXX),icpc)
RELEASE = -O3 -ipo -no-prec-div -xHOST -g
endif
ifeq ($(CXX),g++)
RELEASE = -O3 -g
endif

STDCPP11 = -std=c++11

#CFLAGS = $(DEBUG)
CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)

INCLUDE = -isystem ../../FDPS/src -I../../include

TARGET = mdpd_fluid.out mdpd_fluid_omp.out

all:$(TARGET)

mdpd_fluid.out: mdpd_fluid.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

mdpd_fluid_omp.out: mdpd_fluid.cpp
	$(CXX) $(CFLAGS) $(OPENMP) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
