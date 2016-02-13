CXX = mpicxx

OPENMP = -openmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual -Wunused-variable
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3 -ipo -no-prec-div -xHOST
STDCPP11 = -std=c++11

CFLAGS = $(DEBUG)
#CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)
#

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)
CFLAGS += -DPAIRWISE_DPD
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

INCLUDE = -I../../FDPS/src -I../../

TARGET = io_test_mpi.out

all:$(TARGET)

io_test_mpi.out: io_test.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
