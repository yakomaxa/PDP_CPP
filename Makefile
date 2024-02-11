CXX	= g++-13
CXXFLAGS  =  -std=c++20 -Wall -I./gemmi/include/ -O3 -Wl,-ld_classic
LIBS    = -lm -lz 
VPATH   = ./src/
SRC     = $(shell ls $(VPATH)/*.cpp)
HEADERS = $(shell ls $(VPATH)/*.hpp) 
OBJS    = $(SRC:.cpp=.o)
TARGET  = pdp_2023

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET)
