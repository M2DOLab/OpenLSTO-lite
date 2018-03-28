from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np 
import glob

pyxname = 'pyBind'
sources = glob.glob("./../src/*.cpp")
sources.append(pyxname + '.pyx')

setup(
    ext_modules=cythonize(Extension(
        pyxname, sources=sources,
        language="c++", extra_compile_args=['-std=c++11', '-g','-Wno-sign-compare','-Wno-unused-but-set-variable'],
        include_dirs=[np.get_include(), "../include"],
        library_dirs=['usr/lib', 'usr/', "../include"],
        extra_link_args=["-g"],
    ), gdb_debug=True),

)
