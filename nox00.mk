CXX = icpc
# CXX = mpicxx

ifeq ($(CXX),mpicxx)
use_mpi = yes
endif

# use_omp = yes

MPI = -DPARTICLE_SIMULATOR_MPI_PARALLEL
OPENMP = -DPARTICLE_SIMULATOR_THREAD_PARALLEL
ifeq ($(CXX),icpc)
OPENMP += -openmp
endif
ifeq ($(CXX),mpicxx)
OPENMP += -openmp
endif

WARNINGS = -Wall -Wunused-variable #-Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3

INCLUDE = -I./FDPS/src

# COMMON_FLAGS = $(DEBUG)
COMMON_FLAGS = $(RELEASE)

COMMON_FLAGS += -std=c++11

ifeq ($(CXX),icpc)
OPT_FLAGS = -ipo -no-prec-div
endif
ifeq ($(CXX),mpicxx)
OPT_FLAGS = -ipo -no-prec-div
endif

CXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(OPT_FLAGS)

ifeq ($(use_mpi),yes)
CXX_FLAGS += $(MPI)
endif

ifeq ($(use_omp),yes)
CXX_FLAGS += $(OPENMP)
endif

TARGET_DPD = implicit_dpd_avx.out implicit_dpd_avx2.out implicit_dpd_sse4.out
TARGET_CMAKE = config_maker.out

all:$(TARGET_DPD) $(TARGET_CMAKE)

implicit_dpd_avx.out : main.cpp
	$(CXX) $(CXX_FLAGS) -xAVX $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_avx2.out : main.cpp
	$(CXX) $(CXX_FLAGS) -xCORE-AVX2 $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_sse4.out : main.cpp
	$(CXX) $(CXX_FLAGS) -xSSE4.2 $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_CMAKE): config_maker.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET_DPD) $(TARGET_CMAKE) core.* *~
