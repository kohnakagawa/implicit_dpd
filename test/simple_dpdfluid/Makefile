CXX = g++
#CXX = icpc

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

CFLAGS = $(DEBUG)
#CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
#CFLAGS += $(OPENMP)
CFLAGS += $(STDCPP11)
CFLAGS += -DPAIRWISE_DPD

INCLUDE = -I../../FDPS/src -I../../

OBJECTS = groot_warren.o

TARGET = groot_warren.out

all:$(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(OBJECTS) $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~
