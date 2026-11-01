# This file is for testing purposes and is not meant for production code
CXX      ?= g++
CXXSTD   ?= c++17
CXXFLAGS ?= -std=$(CXXSTD) -O2 -Wall -Wextra
LDFLAGS  ?=

TARGET  := rosenbrock
SRC     := main.cpp
HEADERS := $(wildcard *.hpp)

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) evolution.dat
