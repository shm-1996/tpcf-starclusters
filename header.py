#Data Handling
import numpy as np
import argparse 
import PhysicalConstantsCGS as const    
import pickle
from scipy.optimize import curve_fit

#System modules
import time
from block_timer.timer import Timer
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
mpl.rc_file('/Users/shm/.matplotlib/matplotlibrc',use_default_template=False)


##### Some global quantities #########
arcsec_to_degree = 1./3600.
