CXX = g++
# CXX = icpc
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
ifeq ($(CXX),g++)
OPENMP += -fopenmp
endif
ifeq ($(CXX),mpicxx)
OPENMP += -fopenmp
endif

WARNINGS = -Wall -Wextra -Wunused-variable -Wsign-compare -Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3

INCLUDE = -isystem ./FDPS/src -I./include

# COMMON_FLAGS = $(DEBUG)
COMMON_FLAGS = $(RELEASE)

COMMON_FLAGS += -std=c++11

ifeq ($(CXX),icpc)
OPT_FLAGS = -ipo -no-prec-div -xHOST
endif
ifeq ($(CXX),g++)
OPT_FLAGS = -ffast-math -funroll-loops
endif
ifeq ($(CXX),mpicxx)
OPT_FLAGS = -ffast-math -funroll-loops
endif

CXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(OPT_FLAGS)

ifeq ($(use_mpi),yes)
CXX_FLAGS += $(MPI)
endif

ifeq ($(use_omp),yes)
CXX_FLAGS += $(OPENMP)
endif

TARGET_DPD = implicit_dpd.out
TARGET_DPD_CHEM = implicit_dpd_chem.out
TARGET_CMAKE = config_maker.out
TARGET_ECONF = edit_config.out

all:$(TARGET_DPD) $(TARGET_DPD_CHEM) $(TARGET_CMAKE) $(TARGET_ECONF)

$(TARGET_DPD): ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_DPD_CHEM): ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -DCHEM_MODE $< $(LIBRARY) -o $@

$(TARGET_CMAKE): ./src/config_maker.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_ECONF): ./src/edit_config.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET_DPD) $(TARGET_DPD_CHEM) $(TARGET_CMAKE) $(TARGET_ECONF) core.* *~
