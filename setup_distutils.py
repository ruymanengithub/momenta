#! $HOME/SOFTWARE/anaconda/envs/VISSIM/bin/ python


"""

:History:

2 Nov 2016: created

:author: Ruyman Azzollini (MSSL)
"""

from setuptools import setup, Extension
#from numpy.distutils.core import setup,Extension

import sys

from pdb import set_trace as stop

setup(
    name="momenta",
    version="1.0",
    description="Galaxy SB Moments",
    author="Ruyman Azzollini",
    author_email="ruyman.azzollini@gmail.com",
    url="https://github.com/ruymanengithub/momenta",
    long_description=__doc__,
    packages=['momenta'],
    package_dir={'momenta':'momenta/'},
    package_data={'momenta':['momenta/inputs/*']},
    include_package_data=True,
    zip_safe=False,
)
