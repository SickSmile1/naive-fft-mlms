project('Boussinesq','cpp', version:'1.0.0',default_options:['warning_level=3','cpp_std=c++17'], license: 'MIT')

sourceRoot = meson.source_root()
system = host_machine.system()
architecture = host_machine.cpu_family()
cpp = meson.get_compiler('cpp')

program_sources = ['BoussinesqMain.cpp']

fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: true)
thread_dep = dependency('threads')

project_headers = ['Boussinesq.h',
  'BoussinesqFft.h','BoussinesqMain.h',
  'BoussinesqMlms.h']

project_target = shared_library('')
