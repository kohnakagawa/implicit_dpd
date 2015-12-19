CXX = g++
#CXX = icpc

ifeq ($(CXX),icpc)
OPENMP = -openmp
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp
endif

WARNINGS = -Wall -Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG

ifeq ($(CXX),icpc)
RELEASE = -O3 -xHOST -ipo -no-prec-div
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

INCLUDE = -I./include

OBJECTS = main.o

TARGET = implicit_dpd.out
all:$(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@
.SUFFIXES: .c .o
.c.o:
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< $(LIBRARY) -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(OBJECTS) $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~
