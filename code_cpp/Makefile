# Variables
CXX = g++
CXXFLAGS = -std=c++11
TARGET = main

# WARNING: This code could generate DIFFERENT results when compiled with different CXXFLAGS

# Source files
SRCS = main.cpp utils.cpp constitutive.cpp gradient.cpp config.cpp energy.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to build object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: clean