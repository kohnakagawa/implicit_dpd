#CXX = g++
CXX = icpc

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

#CFLAGS = $(DEBUG)
CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)

INCLUDE = -I../FDPS/src -I../

TARGET = run.out run_omp.out

all:$(TARGET)

run.out: profile_run.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

run_omp.out: profile_run.cpp
	$(CXX) $(CFLAGS) $(OPENMP) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
