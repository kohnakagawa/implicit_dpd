CXX = g++
#CXX = icpc

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
RELEASE = -O3
endif

STDCPP11 = -std=c++11

CFLAGS = $(DEBUG)
#CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)
CFLAGS += -DPAIRWISE_DPD

INCLUDE = -isystem ../../FDPS/src -I../../include

TARGET = bondbend.out

all:$(TARGET)

bondbend.out: bondbendtest.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
