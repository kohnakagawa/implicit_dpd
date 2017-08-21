CXX = g++
MPICXX = mpicxx

OPENMP = -fopenmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual -Wunused-variable
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3 -ipo -no-prec-div -xHOST
STDCPP11 = -std=c++11

CFLAGS = $(DEBUG)
#CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)
CFLAGS += $(OPENMP)
MPI_FLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL

INCLUDE = -isystem../../FDPS/src -I../../include

TARGET = mpi_chem.out chem.out

all:$(TARGET)

mpi_chem.out: chem_test.cpp
	$(MPICXX) -DCHEM_MODE -DLOCAL_CHEM_EVENT $(CFLAGS) $(MPI_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

chem.out: chem_test.cpp
	$(CXX) -DCHEM_MODE -DLOCAL_CHEM_EVENT $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
