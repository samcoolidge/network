# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 22:57:47 2016

@author: samuelcoolidge
"""


# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:15:07 2016

@author: samuelcoolidge
"""

from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import numpy

from distutils.core import setup, Extension 
from Cython.Build import cythonize
import numpy

from Cython.Compiler.Options import directive_defaults

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

extensions = [
    Extension("MC_step_cython" , ["MC_step_cython.pyx"],
              include_dirs=[numpy.get_include()],
              define_macros=[('CYTHON_TRACE', '1')])]

setup(
    ext_modules = cythonize(extensions),
    
)
