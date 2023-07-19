CXX = g++ -O3 -std=c++20  # -Wall -fno-rtti -fno-exceptions -pedantic -lpthread -fno-omit-frame-pointer -DCMAKE_BUILD_TYPE=Release
#CXX = g++ -O3 -std=c++20 -fno-exceptions -fno-rtti
# CXX = -pedantic -fno-elide-constructors -Wall -pedantic -std=c++11
MAIN_BINARIES = $(basename $(wildcard *Main.cpp))
TEST_BINARIES = $(basename $(wildcard *Test.cpp))
HEADERS = $(wildcard *.h) 
OBJECTS = $(addsuffix .o, $(basename $(filter-out %Main.cpp %Test.cpp, $(wildcard *.cpp))))
LIBRARIES = -isystem benchmark/include \
  -Lbenchmark/home/sick/Documents/GIT/benchmark -lbenchmark -lpthread -lfftw3 -lm
# -lfftw3 -lm
# -isystem benchmark\include
.PRECIOUS: %.o
.SUFFIXES:
.PHONY: all compile

all: compile

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

%.o: %.cpp $(HEADERS)
	$(CXX) -c $<

Profile:
	sudo perf record -g ./$(basename $(wildcard *Main.cpp))
	sudo perf report -g 'graph,0.5,caller'
