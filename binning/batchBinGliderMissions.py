import os
import re
import fnmatch
import scipy.io
from binMatFile import binMatFile

gliders = ['henry','grampus']
parentDirs = ['V:/glider/'+g for g in gliders]

missionDirs = []
for pdir in parentDirs:
    allContents = os.listdir(pdir)
    missionDirs += [pdir+'/'+mm+'/' for mm in fnmatch.filter(allContents,'mission[0-9]*')]

for mf in missionDirs:
    mno = re.search('\d*\.*\d',mf).group()
    if re.search('henry',mf):
        glider = 'henry'
    else:
        glider = 'grampus'

    mfile = mf + 'matlabdata/mission'+mno+'.mat'
    outfile = 'V:/Catherine/Glider/DataFiles/fullMissions/'+glider+'-mission'+mno+'.nc'
    if not(os.path.isfile(outfile)):
        try:
            mat = scipy.io.loadmat(mfile)
        except:
            continue
    else:
        continue
        
    fullMission = binMatFile(mat,glider)
    fullMission.to_netcdf(outfile)
