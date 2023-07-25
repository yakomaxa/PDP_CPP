CXX	= g++-12
CXXFLAGS  =  -std=c++17 -Wall -I./  -O0
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
