#CXX = g++
CXX = icpc

WARNINGS = -Wextra -Wunused-variable -Wsign-compare #-Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3
STDCPP11 = -std=c++11

INCLUDE = -I../../FDPS/src -I../../src

ifeq ($(CXX), icpc)
LIBRARY = -L./mdstresslib/build -lmdstress -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
endif
ifeq ($(CXX), g++)
LIBRARY = -L./mdstresslib/build -lmdstress -llapack
endif

# CXX_FLAGS = $(DEBUG)
CXX_FLAGS = $(RELEASE)

CXX_FLAGS += $(STDCPP11)
CXX_FLAGS += $(WARNINGS)

TARGET = post_ana

all:$(TARGET)

$(TARGET): analysis.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~