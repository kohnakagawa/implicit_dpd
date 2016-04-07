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
CFLAGS += -DPAIRWISE_DPD
MPI_FLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL

INCLUDE = -I../../FDPS/src -I../../src

TARGET = mpi_bonded.out bonded.out

all:$(TARGET)

mpi_bonded.out: mpi_bonded.cpp
	$(MPICXX) $(CFLAGS) $(MPI_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

bonded.out: mpi_bonded.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
