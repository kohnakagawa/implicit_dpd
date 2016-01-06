#CXX = g++
CXX = icpc

ifeq ($(CXX),icpc)
OPENMP = -openmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp -DPARTICLE_SIMULATOR_THREAD_PARALLEL
endif

WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual -Wunused-variable
DEBUG = -O0 -g -DDEBUG #-DSANITY_CHECK_REALLOCATABLE_ARRAY

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
CFLAGS += $(OPENMP)
CFLAGS += $(STDCPP11)

INCLUDE = -I./FDPS/src

OBJECTS_DPD = main.o
OBJECTS_CMAKE = config_maker.o

TARGET_DPD = implicit_dpd.out
TARGET_CMAKE = config_maker.out

all:$(TARGET_DPD) $(TARGET_CMAKE)

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@
.SUFFIXES: .c .o
.c.o:
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@

$(TARGET_DPD): $(OBJECTS_DPD)
	$(CXX) $(CFLAGS) $(OBJECTS_DPD) $(LIBRARY) -o $@

$(TARGET_CMAKE): $(OBJECTS_CMAKE)
	$(CXX) $(CFLAGS) $(OBJECTS_CMAKE) $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS_DPD) $(OBJECTS_CMAKE) $(TARGET_DPD) $(TARGET_CMAKE) core.* *~
