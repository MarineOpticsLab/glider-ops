import numpy as np
import pandas as pd
import xarray as xr
import os
import re
import fnmatch
import scipy.io
import datetime as dt
import argparse
from binMatFile import binMatFile

parser = argparse.ArgumentParser(description='''\
  This program bins the glider mission .mat files\
  produced by mdb2matlab.m into netCDF files. The output
  is a file named <glider>-mission<number>.nc 
  e.g. henry-mission1.nc''')

parser.add_argument('--glider', nargs=1, type=str, required=True, choices=['henry','grampus'], help='''\
  Name of the glider, either henry or grampus.''')

parser.add_argument('--mission', nargs=1, type=str, required=True, help='''\
  Number of the mission.''')

args=parser.parse_args()
dict_args=vars(args)
glider = dict_args['glider'][0]
mission = dict_args['mission'][0]
odir = dict_args['odir'][0]

matfile = 'V:/glider/'+glider+'/mission'+mission+'/matlabdata/mission'+mission+'.mat'  #update as needed!

mat = scipy.io.loadmat(matfile,mat_dtype = True)

fullMission = binMatFile(mat,glider)
fullMission.to_netcdf('V:/Catherine/Glider/DataFiles/fullMissions/'+glider+'-mission'+mission+'.nc')  #update as needed!
