CXX = g++ -O3 -std=c++20 -Wall -fno-rtti -fno-exceptions -pedantic -lpthread -fno-omit-frame-pointer
#CXX = g++ -O3 -std=c++20 -fno-exceptions -fno-rtti
# CXX = -pedantic -fno-elide-constructors -Wall -pedantic -std=c++11
MAIN_BINARIES = $(basename $(wildcard *Main.cpp))
TEST_BINARIES = $(basename $(wildcard *Test.cpp))
HEADERS = $(wildcard *.h) 
OBJECTS = $(addsuffix .o, $(basename $(filter-out %Main.cpp %Test.cpp, $(wildcard *.cpp))))
LIBRARIES =  -lbenchmark
# -isystem benchmark\include
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
	rm -f *.gcda
	rm -f perf.data*

%Main: %Main.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBRARIES)

%Test: %Test.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBRARIES) -lgtest -lgtest_main -lpthread

%.o: %.cpp $(HEADERS)
	$(CXX) -c $<

Profile:
	sudo perf record -g ./$(basename $(wildcard *Main.cpp))
	sudo perf report -g 'graph,0.5,caller'
