CXX = g++ -O3 -std=c++20 -mtune=native -march=native -fopenmp# -fsanitize=address
# CXX = g++ -fno-elide-constructors -Wall -pedantic -std=c++11
MAIN_BINARIES = $(basename $(wildcard *Main.cpp))
TEST_BINARIES = $(basename $(wildcard *Test.cpp))
HEADERS = $(wildcard *.h)
OBJECTS = $(addsuffix .o, $(basename $(filter-out %Main.cpp %Test.cpp, $(wildcard *.cpp))))
LIBRARIES = 

.PRECIOUS: %.o
.SUFFIXES:
.PHONY: all compile test

all: compile test

compile: $(MAIN_BINARIES) $(TEST_BINARIES)

test: $(TEST_BINARIES)
	for T in $(TEST_BINARIES); do ./$$T; done


clean:
	rm -f *.o
	rm -f $(MAIN_BINARIES)
	rm -f $(TEST_BINARIES)
	rm -f *.plain

%Main: %Main.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBRARIES)

%Test: %Test.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBRARIES) -lgtest -lgtest_main -lpthread

%.o: %.cpp $(HEADERS)
	$(CXX) -c $<

valgrind:
	valgrind -s --tool=callgrind ./$(basename $(wildcard *Main.cpp))
