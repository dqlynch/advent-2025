CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

all: day11

day11: day11.cpp
	$(CXX) $(CXXFLAGS) -o day11 day11.cpp

run: day11
	time ./day11

clean:
	rm -f day11

.PHONY: all run clean
