from distutils.core import setup, Extension
import numpy
# define the extension module

setup(
    name='twed',
    version='1.0',
    description='Time Warping Edit Distance',
    author='Pierre-François Marteau',
    author_email='pierrefrancois.marteau@gmail.com',
    license='MIT License',
    packages=['twed'],
    install_requires=['numpy'],
    ext_modules=[
        Extension("TWED", ["TWED.c"],
                  include_dirs=[numpy.get_include()]),
    ],
)

