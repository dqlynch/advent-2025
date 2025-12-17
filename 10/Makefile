CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

all: day10

day10: day10.cpp
	$(CXX) $(CXXFLAGS) -o day10 day10.cpp

run: day10
	time ./day10

clean:
	rm -f day10

.PHONY: all run clean
