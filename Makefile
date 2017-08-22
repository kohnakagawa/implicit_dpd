CXX = g++
# CXX = icpc
MPICXX = mpicxx

MPI = -DPARTICLE_SIMULATOR_MPI_PARALLEL
OPENMP = -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

WARNINGS = -Wall -Wextra -Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3

INCLUDE = -isystem ./FDPS/src -I./include

# COMMON_FLAGS = $(DEBUG)
COMMON_FLAGS = $(RELEASE)
COMMON_FLAGS += -std=c++11

ifeq ($(CXX),icpc)
CXX_OPT_FLAGS = -ipo -no-prec-div -xHOST
endif
ifeq ($(CXX),g++)
CXX_OPT_FLAGS = -ffast-math -funroll-loops
endif

MPICXX_OPT_FLAGS = -ffast-math -funroll-loops
# MPICXX_OPT_FLAGS = -ipo -no-prec-div -xHOST

CXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(CXX_OPT_FLAGS)
MPICXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(MPICXX_OPT_FLAGS)

TARGET_DPD_OMP = omp_implicit_dpd.out
TARGET_DPD_OMP_CHEM = omp_implicit_dpd_chem.out
TARGET_DPD_MPI = mpi_implicit_dpd.out
TARGET_DPD_MPI_CHEM = mpi_implicit_dpd_chem.out
TARGET_CMAKE = config_maker.out
TARGET_ECONF = edit_config.out

all:$(TARGET_DPD_OMP) $(TARGET_DPD_OMP_CHEM) $(TARGET_DPD_MPI) $(TARGET_DPD_MPI_CHEM) $(TARGET_CMAKE) $(TARGET_ECONF)

$(TARGET_DPD_OMP): ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(OPENMP) $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_DPD_OMP_CHEM): ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(OPENMP) $(INCLUDE) -DCHEM_MODE $< $(LIBRARY) -o $@

$(TARGET_DPD_MPI): ./src/main.cpp
	$(MPICXX) $(MPICXX_FLAGS) $(MPI) $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_DPD_MPI_CHEM): ./src/main.cpp
	$(MPICXX) $(MPICXX_FLAGS) $(MPI) $(INCLUDE) -DCHEM_MODE $< $(LIBRARY) -o $@

$(TARGET_CMAKE): ./src/config_maker.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_ECONF): ./src/edit_config.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET_DPD_OMP) $(TARGET_DPD_OMP_CHEM) $(TARGET_DPD_MPI_CHEM) $(TARGET_CMAKE) $(TARGET_ECONF) core.* *~
