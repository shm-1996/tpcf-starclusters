#Data Handling
import numpy as np
import argparse 
import PhysicalConstantsCGS as const    
import pickle
from scipy.optimize import curve_fit

#System modules
import time
import os
import sys
import warnings
import shutil
import subprocess

#Astropy Data Handling
from spectral_cube import SpectralCube
from astropy import units as u
import aplpy
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
import scipy

#Visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('classic')
mpl.rc_file('/Users/shm/.matplotlib/matplotlibrc',
    use_default_template=False)


##### Some global quantities #########
arcsec_to_degree = 1./3600.
list_of_galaxies = ['NGC_0628','NGC_1313','NGC_1566','NGC_3344',
        'NGC_3738','NGC_4395','NGC_4449','NGC_4656',
        'NGC_5194','NGC_5474','NGC_6503','NGC_7793']

# Pickle Data Handling
def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(obj, f, protocol=2)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
