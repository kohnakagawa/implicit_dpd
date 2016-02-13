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

INCLUDE = -I../../FDPS/src -I../../

TARGET = scan.out

all:$(TARGET)

scan.out: scan.cpp
	$(MPICXX) $(CFLAGS) $(MPI_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
