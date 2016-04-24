CXX = g++

WARNINGS = -Wextra -Wunused-variable -Wsign-compare #-Wnon-virtual-dtor -Woverloaded-virtual
DEBUG = -O0 -g -DDEBUG
RELEASE = -O3
STDCPP11 = -std=c++11

INCLUDE = -I../FDPS/src -I../src

CXX_FLAGS = $(DEBUG)
#CXX_FLAGS = $(RELEASE)

CXX_FLAGS += $(STDCPP11)
CXX_FLAGS += $(WARNINGS)

TARGET = post_ana

all:$(TARGET)

$(TARGET): analysis.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $< -o $@

clean:
	rm -f $(TARGET) core.* *~