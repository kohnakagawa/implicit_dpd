CXX = icpc
MPICXX = mpicxx

MPI = -DPARTICLE_SIMULATOR_MPI_PARALLEL
OPENMP = -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp

WARNINGS = -Wall -Wunused-variable -Wsign-compare #-Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3

INCLUDE = -I./FDPS/src -I./src

# COMMON_FLAGS = $(DEBUG)
COMMON_FLAGS = $(RELEASE)

COMMON_FLAGS += -std=c++11

OPT_FLAGS = -ipo -no-prec-div

CXX_FLAGS = $(COMMON_FLAGS) $(WARNINGS) $(OPT_FLAGS)

TARGET_DPD = implicit_dpd_avx.out implicit_dpd_avx2.out implicit_dpd_sse4.out implicit_dpd_avx_mpi.out implicit_dpd_avx2_mpi.out implicit_dpd_sse4_mpi.out
TARGET_CMAKE = config_maker.out

all : $(TARGET_DPD) $(TARGET_CMAKE)

implicit_dpd_avx2.out : ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(OPENMP) -xCORE-AVX2 $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_avx.out : ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(OPENMP) -xAVX $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_sse4.out : ./src/main.cpp
	$(CXX) $(CXX_FLAGS) $(OPENMP) -xSSE4.2 $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_avx2_mpi.out : ./src/main.cpp
	$(MPICXX) $(CXX_FLAGS) $(MPI) -xCORE-AVX2 $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_avx_mpi.out : ./src/main.cpp
	$(MPICXX) $(CXX_FLAGS) $(MPI) -xAVX $(INCLUDE) $< $(LIBRARY) -o $@

implicit_dpd_sse4_mpi.out : ./src/main.cpp
	$(MPICXX) $(CXX_FLAGS) $(MPI) -xSSE4.2 $(INCLUDE) $< $(LIBRARY) -o $@

$(TARGET_CMAKE): ./src/config_maker.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET_DPD) $(TARGET_CMAKE) core.* *~
