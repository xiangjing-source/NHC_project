CXX = g++
CXXFLAGS = -O2 -I./include
SRC = src/main.cpp src/io.cpp src/force.cpp src/thermostat_nhc.cpp src/globals.cpp
OBJ = $(patsubst src/%.cpp,bin/%.o,$(SRC))
TARGET = md

all: bin $(TARGET)

bin:
	mkdir -p bin

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Compilation successful!"  
bin/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf bin $(TARGET)
