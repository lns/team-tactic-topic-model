VERSION = 1

CXX = g++
CXXFLAGS = -Wall -std=c++0x -O2
#CXXFLAGS = -Wall -std=c++11 -O0 -Og -g


all: main

main: main.cpp lda.cpp rng.cpp pct.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) main

