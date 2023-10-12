# Boussinesq Kernel Calculations with Naive, Multi-Level-Multi-Summation and FFT implementation #

# Requirements #

<p>FFTW3 <br>
Google Test Framework - https://github.com/google/googletest <br>
Google Benchmark Framework - https://github.com/google/benchmark <br>
</p>

<p>Note: Google Test and Benchmark Framework are available from: <br>
<ul>
  <li>apt: sudo apt install libgtest-dev libbenchmark-dev</li>
  <li>brew: brew install googletest google-benchmark</li>
</ul>
</p>

# Makefile #

## Options ## 
<p>
make all - compiles all and calls Tests and Valgrind with cachegrind (further installations are necessary - https://valgrind.org/)
</p>

By default, a call to <b>make</b> in the root directory of the repo is sufficient. If Google Test and Benchmark are not installed, the call can be modified to <b>make compile</b> to skip the tests.

# Bench #

# Tests #
