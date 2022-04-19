import xarray as xr
import numpy as np
import pandas as pd
import os 
import re


parentDir = '../../data/mission-netcdf/' #update as needed!
files = os.listdir(parentDir)
gFiles = [parentDir + ff for ff in files if os.path.isfile(parentDir+ff)]

for gf in gFiles:
    gdata = xr.open_dataset(gf)
    deepchl1 = gdata['chlfl_east'].sel(zbin=slice(150,180,1))
    deepchl2 = gdata['chlfl_west'].sel(zbin=slice(150,180,1))
    deepchl = np.concatenate((deepchl1.data,deepchl2.data))
    offset = np.nanmedian(deepchl)
    
    for dd in ['east','west']:
        gdata['chlfl_'+dd] = gdata['chlfl_'+dd] - offset
    
    newfile = '../../data/mission-netcdf/darkCountCorrected/'+re.split('/',gf)[-1]  #update as needed!
    gdata.to_netcdf(newfile)