# Boussinesq Kernel Calculations with Naive, Multi-Level-Multi-Summation and FFT implementation #

# Requirements #

<p>FFTW3 <br>
Google Test Framework - https://github.com/google/googletest <br>
Google Benchmark Framework - https://github.com/google/benchmark <br>
python3 <br>
</p>

<p>Note: Google Test and Benchmark Framework are available from: <br>
<ul>
  <li>apt: sudo apt install libgtest-dev libbenchmark-dev</li>
  <li>brew: brew install googletest google-benchmark</li>
</ul>
</p>

# Makefile #

By default, a call to <b>make</b> in the root directory of the repo is sufficient. If Google Test and Benchmark are not installed, the call can be modified to <b>make compile</b> to skip the tests.

## Options ## 
<p>
  <ul>
    <li>make all - compiles all and calls Tests and Valgrind with cachegrind (further installations are necessary - https://valgrind.org/)</li>
    <li>make compile - skip tests and valgrind </li>
    <li>make test - run tests only</li>
    <li>make valgrind - run memory leak checks and see assembly, lots of options to choose from: https://valgrind.org/</li>
    <li>make clean - clean *.o files and other residue from compilation </li>
    <li>make checkstyle - run static checker cpplint for project</li>
    <li>make bench - compile benchmark</li>
  </ul>
</p>

# Source of Algorithms #
https://www.sciencedirect.com/science/article/abs/pii/S0043164899001131
https://iopscience.iop.org/article/10.1088/2051-672X/ac860a/pdf
