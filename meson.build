project('Boussinesq','cpp', 
  version:'1.0.0',
  default_options:['warning_level=3','cpp_std=c++17','buildtype=debugoptimized'],
  license: 'MIT'
  )

project_headers = []

sourceRoot = meson.source_root()
system = host_machine.system()
architecture = host_machine.cpu_family()
cpp = meson.get_compiler('cpp')

inc = subdir('inc')

gtest = subproject('gtest')
benchmarks = subproject('google-benchmark')

fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: true)
thread_dep = dependency('threads')
gtest_dep = gtest.get_variable('gtest_dep') 
gbenchmark_dep = benchmarks.get_variable('google_benchmark_dep')


sources = ['BoussinesqMain.cpp','BoussinesqMlms.cpp', 'BoussinesqFft.cpp', 'Boussinesq.cpp']
my_executable = executable('BoussinesqMain', sources, include_directories:inc, dependencies:[fftw3_dep])

sources2 = ['BoussinesqTest.cpp','BoussinesqMlms.cpp', 'BoussinesqFft.cpp', 'Boussinesq.cpp']
my_executable2 = executable('BoussinesqTest', sources2, include_directories:inc, dependencies:[fftw3_dep, thread_dep,gtest_dep])

sources3 = ['BenchMain.cpp','BoussinesqMlms.cpp', 'BoussinesqFft.cpp', 'Boussinesq.cpp']
my_executable3 = executable('Bench', sources3, include_directories:inc, dependencies:[fftw3_dep, thread_dep,gbenchmark_dep])
# project_target = shared_library('')
