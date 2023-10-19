project('Boussinesq','cpp', 
  version:'1.0.0',
  default_options:['warning_level=3','cpp_std=c++17','buildtype=debugoptimized'],
  license: 'MIT'
  )

cpplint_program = find_program('cpplint', required : false)

if cpplint_program.found()
    # - https://github.com/cpplint/cpplint

    cpplint_config = ['--verbose=0',
      '--linelength=130',
      '--filter=-build/include_subdir,-build/header_guard,-legal/copyright,'+
      '-build/include_what_you_use,-readability/casting,-build/include_subdir,'+
      '-build/include_order,-readability/multiline_string',
      '--verbose=5',
      '--recursive',
      join_paths(meson.current_source_dir(), 'src')]
    cpplint_target = run_target('cpplint', command : [cpplint_program, cpplint_config])
endif

sourceRoot = meson.source_root()
system = host_machine.system()
architecture = host_machine.cpu_family()
cpp = meson.get_compiler('cpp')

inc = include_directories('headers')

gtest = subproject('gtest')
benchmarks = subproject('google-benchmark')

fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: true)
thread_dep = dependency('threads')
gtest_dep = gtest.get_variable('gtest_dep')
gtest_main_dep = gtest.get_variable('gtest_main_dep')
gbenchmark_dep = benchmarks.get_variable('google_benchmark_dep')


sources = ['BoussinesqMain.cpp', 'src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable = executable('BoussinesqMain', sources, include_directories:inc, dependencies:[fftw3_dep])

sources2 = ['BoussinesqTest.cpp','src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable2 = executable('BoussinesqTest', sources2, include_directories:inc, dependencies:[fftw3_dep, thread_dep,gtest_dep,gtest_main_dep])

sources3 = ['BenchMain.cpp','src/BoussinesqMlms.cpp', 'src/BoussinesqFft.cpp', 'src/Boussinesq.cpp']
my_executable3 = executable('Bench', sources3, include_directories:inc, dependencies:[fftw3_dep, thread_dep,gbenchmark_dep])
# project_target = shared_library('')
