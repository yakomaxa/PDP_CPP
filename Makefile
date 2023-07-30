CXX	= g++-12
CXXFLAGS  =  -std=c++20 -Wall -I./gemmi/include/ -O3 -ffast-math -Ofast -march=native -funroll-loops -flto
LIBS    = -lm -lz 
VPATH   = ./
SRC     = $(shell ls $(VPATH)/*.cpp)
HEADERS = $(shell ls $(VPATH)/*.hpp) 
OBJS    = $(SRC:.cpp=.o)
TARGET  = pdp_2023

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET)
