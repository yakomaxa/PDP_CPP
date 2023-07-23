#CXX	= /lab/home/ksakuma/work/build_gcc11/gcc-11.3.0/bin/g++-11.3.0
CXX	= g++-12
CXXFLAGS  =  -std=c++20 -Wall -I./ 
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
