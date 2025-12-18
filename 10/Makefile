CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

all: part1 part2

part1: part1.cpp
	$(CXX) $(CXXFLAGS) -o part1 part1.cpp

part2: part2.cpp
	$(CXX) $(CXXFLAGS) -o part2 part2.cpp

run1: part1
	time ./part1

run2: part2
	time ./part2

clean:
	rm -f part1 part2

.PHONY: all run clean
