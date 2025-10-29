from distutils.core import setup, Extension
import numpy
# define the extension module
TWED = Extension(
    'TWED', sources=['TWED.c'],
    include_dirs=[numpy.get_include()],
    )

# run the setup
setup(ext_modules=[TWED])


