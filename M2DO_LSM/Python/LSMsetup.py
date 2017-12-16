from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np 
import glob

pyxname = 'py_lsmBind'
sources = glob.glob("./../src/*.cpp")
sources.remove("./../src/sensitivity.cpp")
sources.remove("./../src/optimise.cpp")
sources.append(pyxname + '.pyx')

setup(
    ext_modules=cythonize(Extension(
        pyxname, sources=sources,
        language="c++", extra_compile_args=['-std=c++14','-D PYBIND', '-g','-Wno-sign-compare','-Wno-unused-but-set-variable'],
        include_dirs=[np.get_include(), "../include"],
        library_dirs=['usr/lib', 'usr/', "../include"],
        extra_link_args=["-g"],
    ), gdb_debug=True),

)
