from distutils.core import setup, Extension
import numpy
# define the extension module
TWED = Extension('TWED', sources=['TWED.c'])

# run the setup
setup(ext_modules=[TWED])

