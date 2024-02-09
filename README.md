# naive-fft-mlms
mlms - rough surfaces
## add pybind11
git submodule add -b stable ../../pybind/pybind11 extern/pybind11
git submodule update --init
## add pretty-print for eigen
git clone https://github.com/dmillard/eigengdb
paste this line inside bin/eigengdb_register_printers
  loader2 = 'python\nimport numpy as np\nnp.set_printoptions(linewidth=200, formatter={"float": lambda x: "{:5.2f}".format(x) if x !=0 else "     "})\nend'
and add   f.write(loader2)
python3 -m pip install eigengdb # Make sure to use system python which matches GDB
python3 bin/eigengdb_register_printers
