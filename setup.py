from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('gl_cont_cython.pyx'))