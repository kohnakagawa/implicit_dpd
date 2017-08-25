CXX = g++

CFLAGS = -std=c++11 -g -O3 -ffast-math -funroll-loops -pthread
# CFLAGS = -std=c++11 -O0 -g

OBJECTS = sysdraw.o mousehandle.o jpegout.o main.o
TARGET = cDraw
LOADLIBES = -lpthread -lglut -lGLU -lGL -ljpeg -lGLEW

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(OBJECTS) $(LOADLIBES) -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~
