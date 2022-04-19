import xarray as xr
import numpy as np
import os
import re

def signChangeLoc(z):
    ''' returns the 1st index of the input array where 
    the elements in the input array change sign
    input: z = a profile'''
    
    nonans = ~np.isnan(np.sign(z))
    signChanges = np.where(np.abs(np.diff(np.sign(z[nonans])))>1)[0]
    if signChanges.size == 0:
        signChange = np.nan
    else:
        znonans = np.where(nonans)
        signChange = float(znonans[0][signChanges[0]+1])
    return(signChange)

def findingNights(gdata):
    '''finding nights based on surface Ed490
    inputs: gdata = xarray dataset containing all the glider data'''
    
    surfdata = gdata.where(gdata.zbin<10).mean(dim='zbin')
    nightval = 0.1

    nightData = {}
    for direction in ['east','west']:
        night_index = surfdata['Ed490_'+direction] < nightval
        nights = surfdata['Ed490_'+direction][night_index]

        if len(nights) > 0:
            nData = gdata.where(nights)

            nightData[direction] = nData[[vv for vv in gdata.variables if re.search('_'+direction,vv)]]

        else:
            nightData[direction] = gdata[[vv for vv in gdata.variables if re.search('_'+direction,vv)]]*np.nan        
    
    return nightData

def refArrayMeanNights(nightData,gdata,direction,lons):
    '''finding the mean chlfl, bb and chlfl:bb for a full night and 
    interpolating between each night
    inputs: nightData = xarray dataset containing only the night profiles (output of findingNights)
            direction = string, either 'east' or 'west'
            lons = the longitudes to interpolate between each night (= gdata.lonbins)'''
    
    #finding the edges of each night (i.e. where night turns to day)
    #based on longitude differences
    nightLonDiff = np.diff(nightData.lonbin)
    nidx = np.where(nightLonDiff>0.1)[0]+1
    nidx = np.concatenate(([0],nidx))
    sLon = nightData.lonbin[nidx].data
    eLon = np.concatenate((nightData.lonbin[nidx[1:]-1],[nightData.lonbin[-1]]))

    #setting some parameters
    variables = ['chlfl','bb532']

    nightArrays = {}
    for vv in variables:
        varndata = nightData[vv+'_'+direction]

        #calculating the average chl profile for each night
        nightav = []
        for ii in range(0,len(nidx)-1):
            nightav += [varndata[:,nidx[ii]:nidx[ii+1]].mean(dim='lonbin')]
        nightav += [varndata[:,nidx[-1]:].mean(dim='lonbin')]

        #making a reference DataArray of the night means assigned to the 297 lonbins
        varnight = []
        for il in lons:
            dix = True 
            for ix,(ils,ile) in enumerate(zip(sLon,eLon)):
                if (il >= ils) & (il <= ile): 
                    varnight += [nightav[ix]] 
                    dix = False 
                    break 
            if dix: 
                varnight += [nightav[ix]*np.nan]

        # converting the list to a DataArray
        varNightArray = xr.concat(varnight,dim='lonbin')
        varNightArray.coords['lonbin'] = lons

        # interpolating the night means over the day segments
        varRefArray = varNightArray.interpolate_na(dim='lonbin')

        nightArrays[vv+"_"+direction] = varRefArray.T
    refArrays = xr.Dataset(nightArrays)
    refArrays['chlbb_'+direction] = refArrays['chlfl_'+direction]/refArrays['bb532_'+direction]
    
    return refArrays

def refArrayEdgeNights(nightData,gdata,direction,lons):
    '''finding the first/last chlfl, bb and chlfl:bb for a night and 
    interpolating between each night (over the day)
    inputs: nightData = xarray dataset containing only the night profiles (output of findingNights)
            gdata = xarray dataset containing night and day profiles
            direction = string, either 'east' or 'west'
            lons = the longitudes to interpolate between each night (= gdata.lonbins)'''
    
    #finding the edges of each night (i.e. where night turns to day)
    #based on longitude differences
    nightLonDiff = np.diff(nightData.lonbin)
    nidx = np.where(nightLonDiff>0.1)[0]+1
    #nidx = np.concatenate(([0],nidx))
    sLon = nightData.lonbin[nidx].data
    eLon = nightData.lonbin[nidx-1].data

    #setting some parameters
    variables = ['chlfl','bb532']

    nightArrays = {}
    for vv in variables:
        varndata = nightData[vv+'_'+direction]

        #making a reference DataArray of the first/last nights assigned to the 297 lonbins
        varnight = gdata[vv+'_'+direction].data.copy()
        for ix,il in enumerate(lons):
            for iy,(ils,ile) in enumerate(zip(sLon,eLon)):
                if (il <= ils) & (il >= ile): 
                    varnight[:,ix] = varndata[:,0]*np.nan
                    break

        # converting the list to a DataArray
        varNightArray = xr.DataArray(varnight,coords = [gdata.zbin,lons], dims=['zbin','lonbin'])

        # interpolating the night means over the day segments
        varRefArray = varNightArray.interpolate_na(dim='lonbin')

        nightArrays[vv+"_"+direction] = varRefArray
    refArrays = xr.Dataset(nightArrays)
    refArrays['chlbb_'+direction] = refArrays['chlfl_'+direction]/refArrays['bb532_'+direction]

    return refArrays

def findQuenchingDepthSC(gdata,refArrays,direction,zeu=None):
    '''quenching depth using sign change'''
    chldiff = gdata['chlfl_'+direction] - refArrays['chlfl_'+direction]
    qdepths = np.apply_along_axis(signChangeLoc, 0, chldiff)
    qdepths[qdepths < 5] = np.nan
    
    return qdepths

def findQuenchingDepthSCInEZ(gdata,refArrays,direction,zeu):
    '''quenching depth using sign change
    in the euphotic zone'''
    chldiff = gdata['chlfl_'+direction] - refArrays['chlfl_'+direction]
    qdepths = np.apply_along_axis(signChangeLoc, 0, chldiff)
    qdepths[qdepths < 5] = np.nan
    qdepths[qdepths > zeu] = np.nan
       
    return qdepths

def findQuenchingDepthSCupper(gdata,refArrays,direction,zeu=None):
    '''quenching depth using sign change
    in the uppper 65m'''
    chldiff = gdata['chlfl_'+direction] - refArrays['chlfl_'+direction]
    qdepths = np.apply_along_axis(signChangeLoc, 0, chldiff)
    qdepths[qdepths < 5] = np.nan
    qdepths[qdepths > 65] = np.nan
       
    return qdepths

def correctQuenching(gdata,refArrays,direction,qdepths):
    #correcting quenching
    surfdata = gdata.where(gdata.zbin<10).mean(dim='zbin')
    nightval = 0.1
    night_index = surfdata['Ed490_'+direction] < nightval
    
    days = surfdata['Ed490_'+direction][~night_index]
    dayData = gdata.where(days)
    refData = refArrays.where(days)
    qdep = qdepths[~night_index]
    chlCorr = []
    for ii in range(0, len(days)):
        try:
            qz = int(qdep[ii])
            chl = dayData['chlfl_'+direction][:,ii]
            bb = dayData['bb532_'+direction][:,ii]
            nchlbb = refData['chlbb_'+direction][:,ii]

            #initial correction
            inicorr = chl.where(chl.zbin>qz+1,nchlbb*bb)

            #2nd correction
            chlCorr += [inicorr.where((inicorr-chl)>0,chl)]
        except ValueError:
            chlCorr += [dayData['chlfl_'+direction][:,ii]]
    chlCorrArray = xr.concat(chlCorr,dim='lonbin')

    newChl = chlCorrArray.combine_first(gdata['chlfl_'+direction])
    
    return newChl.T

def findingDays(gdata):
    '''finding days based on surface Ed490
    This function is a copy of findingNights, but with the Ed490
    condition switched around
    inputs: gdata = xarray dataset containing all the glider data'''
    
    surfdata = gdata.where(gdata.zbin<10).mean(dim='zbin')
    nightval = 0.1

    dayData = {}
    for direction in ['east','west']:
        day_index = surfdata['Ed490_'+direction] > nightval
        days = surfdata['Ed490_'+direction][night_index]

        if len(days) > 0:
            nData = gdata.where(days)

            dayData[direction] = nData[[vv for vv in gdata.variables if re.search('_'+direction,vv)]]

        else:
            dayData[direction] = gdata[[vv for vv in gdata.variables if re.search('_'+direction,vv)]]*np.nan        
    
    return dayData


    