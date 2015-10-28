from distutils.core import setup, Extension
import numpy.distutils.misc_util
 
module1 = Extension ('tv', 
                     library_dirs = ['/usr/lib/'],
                     include_dirs = ['/usr/include/python2.6/','/usr/share/pyshared/numpy/core/include','/usr/include/', numpy.distutils.misc_util.get_numpy_include_dirs()],
                     libraries = ['gsl', 'gslcblas'],
                     sources = ['tvmodule.c', 'tvmethods.c']
                     )

setup (name = 'tv',
       version = '1.0',
       description = 'TV package',
       ext_modules = [module1]
       )
