from distutils.core import setup
from Cython.Build import cythonize
import numpy

# python density_help_setup.py build_ext --inplace

setup(
    ext_modules = cythonize("*.pyx", 
                            annotate = False,
                            language_level = "3"),
    include_dirs=[numpy.get_include()]
)



