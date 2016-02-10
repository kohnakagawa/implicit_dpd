CXX = g++
#CXX = icpc

#use_tbb = yes

ifeq ($(CXX),icpc)
OPENMP = -openmp
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp
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

#CFLAGS = $(DEBUG)
CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)

CFLAGS += $(WARNINGS)
CFLAGS += $(STDCPP11)
CFLAGS += $(OPENMP)

ifeq ($(use_tbb),yes)
LIBRARY = -ltbb
endif

TARGET = sort.out

all:$(TARGET)

sort.out: main.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
