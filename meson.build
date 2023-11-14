project('Boussinesq','cpp', 
  version:'1.0.0',
  default_options:['warning_level=3','cpp_std=c++20','buildtype=debug'],
  license: 'MIT'
  )


##
# cpp-lint
##
cpplint_program = find_program('cpplint.py', required : false)
if cpplint_program.found()
    # - https://github.com/cpplint/cpplint

    cpplint_config = ['--verbose=1',
      '--linelength=130',
      '--filter=-build/include_subdir,-build/header_guard,-legal/copyright,'+
      '-build/include_what_you_use,-readability/casting,-build/include_subdir,'+
      '-build/include_order,-readability/multiline_string',
      '--verbose=5',
      '--recursive',
      join_paths(meson.current_source_dir(), 'src')]
    cpplint_target = run_target('cpplint.py', command : [cpplint_program, cpplint_config])
endif

##
# python module and pybind11 dep
##
py3_inst = import('python').find_installation('python3', required:true)
py11_dep = include_directories('extern/pybind11/include')

sourceRoot = meson.source_root()
system = host_machine.system()
architecture = host_machine.cpu_family()
cpp = meson.get_compiler('cpp')

inc = include_directories('headers')

##
# Gtest and GBench framworks, Eigen from wrapDB, 
##
gtest = subproject('gtest')
benchmarks = subproject('google-benchmark')
eigen = subproject('eigen')

# get dependencies
fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: true)
thread_dep = dependency('threads')
eigen_dep = eigen.get_variable('eigen_dep')
gtest_dep = gtest.get_variable('gtest_dep')
gtest_main_dep = gtest.get_variable('gtest_main_dep')
gbenchmark_dep = benchmarks.get_variable('google_benchmark_dep')
omp = dependency('openmp')

error_args = ['-Wall', '-Werror', '-Wshadow', '-Wundef','-fno-common','-g3', '-O0','-pedantic','-fno-omit-frame-pointer']

# compile main
sources = ['BoussinesqMain.cpp', 'src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable = executable('BoussinesqMain', sources, include_directories:inc, dependencies:[fftw3_dep, eigen_dep,omp])#, cpp_args: error_args)

# compile Tests
sources2 = ['BoussinesqTest.cpp','src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable2 = executable('BoussinesqTest', sources2, include_directories:inc,
  dependencies:[fftw3_dep, thread_dep,gtest_dep,gtest_main_dep, eigen_dep,omp], cpp_args: error_args)
test('gtest',my_executable2)

# compile Benchmark
sources3 = ['BenchMain.cpp','src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable3 = executable('Bench', sources3, include_directories:inc, dependencies:[fftw3_dep, thread_dep,gbenchmark_dep, eigen_dep,omp], cpp_args: error_args)

#my_exe4 = executable('sparseTest',['testfolder/test.cpp'], dependencies:[thread_dep, gbenchmark_dep, eigen_dep])
#my_exe5 = executable('ccg',['testfolder/Ccg.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp'], include_directories:inc, dependencies:[eigen_dep, fftw3_dep])

# create python package for Mlms, Fft
py3_inst.extension_module('boussinesq',
  sources: ['BoussinesqBind.cpp'],
  include_directories: [py11_dep, inc],
  dependencies: [py3_inst.dependency(), eigen_dep],
)
