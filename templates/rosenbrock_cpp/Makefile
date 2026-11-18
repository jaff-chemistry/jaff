# Builds the example driver (main.cpp) against the header-only integrator.
# A host code only needs to include solver.hpp (it pulls in chemistry.hpp,
# user_functions.hpp and tableau.hpp); no library is built.
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall
TARGET   := test
HEADERS  := chemistry.hpp user_functions.hpp solver.hpp tableau.hpp

.PHONY: all run clean

all: $(TARGET)

$(TARGET): main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ main.cpp

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) solution.dat
