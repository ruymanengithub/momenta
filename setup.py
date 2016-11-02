#! $HOME/SOFTWARE/anaconda/envs/VISSIM/bin/ python


"""

:History:

2 Nov 2016: created

:author: Ruyman Azzollini (MSSL)
"""

from numpy.distutils.core import setup,Extension
from numpy.distutils.misc_util import Configuration

import sys
import os

from pdb import set_trace as stop


def configuration():
    
    config = Configuration()
    
    config.add_subpackage('momenta')
    
    config.add_data_dir(['momenta/inputs','momenta/inputs'])
    config.add_data_dir(['momenta/doc','momenta/doc'])
       
    return config


def setup_package():
        
    metadata = dict(
     name="momenta",
     version="1.0",
     description="Galaxy SB Momenta",
     author="Ruyman Azzollini",
     author_email="ruyman.azzollini@gmail.com",
     url="https://github.com/ruymanentighub/momenta",
     configuration=configuration)
    
    setup(**metadata)
    
    return
    

if __name__ == '__main__': setup_package()

