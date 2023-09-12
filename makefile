CXX = g++ -g -O3 -std=c++17 -fno-rtti -fno-exceptions -pedantic -lpthread #-fsanitize=address-DNDEBUG #-mtune=native -march=native -fopenmp# 
# CXX = g++ -fno-elide-constructors -Wall -pedantic -std=c++11
MAIN_BINARIES = $(basename $(wildcard *Main.cpp))
TEST_BINARIES = $(basename $(wildcard *Test.cpp))
HEADERS = $(wildcard *.h)
OBJECTS = $(addsuffix .o, $(basename $(filter-out %Main.cpp %Test.cpp, $(wildcard *.cpp))))
LIBRARIES = -lfftw3 -lm 

.PRECIOUS: %.o
.SUFFIXES:
.PHONY: all compile test checkstyke

all: compile test checkstyle

compile: $(MAIN_BINARIES) $(TEST_BINARIES)

test: $(TEST_BINARIES)
	for T in $(TEST_BINARIES); do ./$$T; done

checkstyle:
	python3 cpplint.py --repository=. *.h *.cpp

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
	valgrind -s --tool=callgrind --dump-instr=yes --collect-jumps=yes ./$(basename $(wildcard *Main.cpp))
