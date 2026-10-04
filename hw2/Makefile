CXX = g++ # Replace with your preferred compiler or keep g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++17 
# I wanted to do c++20 but VSC somehow doesn't let me and I dont want to debug with the makefile so c++17 it is.
# I am using 03 to optimize for execution speed.
# I am also using Wall and Wextra to see if I have any issues with the code

.PHONY: all carbon golf clean

all: carbon golf

# Compiles the binary if needed using the above compiler and flags, then executes with or without WIDTH
carbon: carbon_bin
	./carbon $(if $(WIDTH),--plot=$(WIDTH),)

carbon_bin: carbon.cpp
	$(CXX) $(CXXFLAGS) -o carbon carbon.cpp

# Compiles the binary if needed, then executes with or without THETA
golf: golf_bin
	./golf $(if $(THETA),--plot=$(THETA),)

golf_bin: golf.cpp
	$(CXX) $(CXXFLAGS) -o golf golf.cpp

clean:
	rm -f carbon golf *.o