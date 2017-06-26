import sys
import os
import re
from setuptools import setup
from setuptools import Extension
#from Cython.Build import cythonize
import numpy
import glob
from setuptools import find_packages
'''
'biopython',\
        'pymc>=2.3.4, < 3.0.0',\
        'scikit-learn>=0.15.2, <= 0.16.1',\
        'statsmodels>=0.5.0',\
        'mpmath>=0.19',\
        'pandas>=0.16.0',\
        'argparse',\
        'numpy',\
'''

input_data_list_commands = glob.glob('mpathic_tests/commands/*.txt')
input_data_list_inputs = glob.glob('mpathic_tests/input/*')

# DON'T FORGET THIS
ext_modules = Extension("mpathic_for_cluster.fast",["src/fast.c"])

# main setup command
setup(
    name = 'mpathic', 
    description = 'Tools for analysis of Sort-Seq experiments.',
    version = '0.01.13',
    author = 'Bill Ireland',
    author_email = 'wireland@caltech.edu',
    #long_description = readme,
    install_requires = [\
        'scipy'
        ],
    platforms = 'Linux (and maybe also Mac OS X).',
    packages = ['mpathic_for_cluster'] + find_packages(),
    package_dir = {'mpathic_for_cluster':'src'},
    scripts = [
            'scripts/mpathic_for_cluster'
            ],
    zip_safe=False,
    ext_modules = [ext_modules],
    include_dirs=['.',numpy.get_include()],
    #package_data = {'mpathic':['tests/*']} # data for command line testing
)


