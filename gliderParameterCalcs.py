import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import savgol_filter
import re

def dfdx(f,x):
    """
    Calculating the derivative of f w.r.t. x
    """
    #making sure the variables are floats
    f = f.astype(float)
    x = x.astype(float)
    
    fprime = np.ones(f.shape)*np.nan
    
    for ii in range(1,len(f)-1):
        #for the polynomial interpolation method
        xdi = x[ii+1] - x[ii]
        xdi_1 = x[ii] - x[ii-1]
        qdx1 = -xdi / (xdi_1*(xdi_1 + xdi))
        qdx2 = (xdi - xdi_1) / (xdi_1 * xdi)
        qdx3 = xdi_1 / (xdi*(xdi_1 + xdi))

        fprime[ii] = qdx1*f[ii-1] + qdx2*f[ii] + qdx3*f[ii+1]
    fprime[0] = (f[1]-f[0]) / (x[1]-x[0])
    fprime[-1] = (f[-1]-f[-2]) / (x[-1]-x[-2])
    
    return fprime

def smoothing(gdata,profilevar,window,polyorder):
    """
    interpolates and smooths glider profiles

    gdata = raw data (xarray dataset)
    profilevar = name of variable in gdata to smooth (string)
    window = window length for smoothing (int)
    polyorder = order of polynomial for smoothing (int)
    """

    interpolated = gdata.interpolate_na(dim='zbin')

    var_smooth = {}
    directions = ['east','west']
    for dd in directions:
        vardata = interpolated[profilevar+'_'+dd].data
        depth = interpolated.zbin

        var1d = np.ones(vardata.shape)*np.nan
        for ix in range(0,297):
            vd = vardata[:,ix]
            nans = np.isnan(vd)
            if len(vd[~nans]) == 0:
                continue
            else:
                try:
                    vdsmooth = savgol_filter(vd[~nans],window,polyorder)
                    vd[~nans] = vdsmooth
                except ValueError:
                    vd[~nans] = np.ones(np.sum(~nans))*np.nan
                var1d[:,ix] = vd

        var_smooth[dd] = var1d.squeeze()

    return var_smooth

def surfDepthChange(chl_smooth,bb_smooth):
    chlDiff = {}
    bbDiff = {}
    zpeak = {}
    depth = np.arange(1,201)
    for dd in ['east','west']:
        chls = chl_smooth[dd]
        bbs = bb_smooth[dd]
        #identifying peaks + calculating %change
        z = np.ones((1,297))*np.nan
        chlChange = np.ones((1,297))*np.nan
        bbChange = np.ones((1,297))*np.nan
        for ip in range(0,297):
            chlProfile = chls[:,ip]
            bbProfile = bbs[:,ip]
            nans = np.isnan(chlProfile)
            chlProfilenonan = chlProfile[~nans]
            nans = np.isnan(bbProfile)
            bbProfilenonan = bbProfile[~nans]
            try:
                chlpeak = np.nanargmax(chlProfile)
                if len(chlProfilenonan) == 0:
                    chlChange[0,ip] = np.nan
                else:
                    surf = chlProfilenonan[0]
                    chlChange[0,ip] = ((chlProfile[chlpeak]-surf)/chlProfile[chlpeak])*100
                    if (chlChange[0,ip] > 10) & (depth[chlpeak] < 80):
                        z[0,ip] = depth[chlpeak]
                if len(bbProfilenonan) == 0:
                    bbChange[0,ip] = np.nan
                else:
                    surf = bbProfilenonan[0]
                    bbChange[0,ip] = ((bbProfile[chlpeak]-surf)/bbProfile[chlpeak])*100
            except ValueError:
                z[0,ip] = np.nan
                chlChange[0,ip] = np.nan
                bbChange[0,ip] = np.nan

        chlDiff[dd] = chlChange
        bbDiff[dd] = bbChange
        zpeak[dd] = z
    
    return chlDiff, bbDiff, zpeak

def PARcalc(gdata):
    #constants for converting to photon number
    h = 6.63 * 10**-34
    c = 3 * 10**8
    NA = 6.02 * 10**23

    PAR = {}
    for dd in ['east','west']:
        #integration (using right hand edge of band for Ed value)
        regex = re.compile('Ed\d\d\d_'+dd)
        Edvars = []
        Edlist = []
        for vv in gdata.variables:
            m = regex.search(vv)
            if m:
                Edvars += [vv]
                Edlist += [gdata[vv].data]
            
        Ed = np.asarray(Edlist[1:])/100

        #manipulating wavelengths for integration
        if len(Edvars)<6:
            wl = np.array([380,443,490,532])
        else:
            wl = np.array([412,443,490,510,532,555,670])
        dwl = np.diff(wl)
        dwl2 = dwl/2
        bandCentre = np.asarray([w + dw for w, dw in zip(wl[:-1],dwl2)])

        par = np.ones((200,297))*np.nan
        Ed = np.swapaxes(Ed.T,0,1)
        for ir,row in enumerate(Ed):
            for ic,col in enumerate(row):
                par[ir,ic] = np.asarray([np.sum((col*dwl)*(bandCentre*10**-9/(h*c)))])
        
        par = (np.asarray(par)/NA)*10**6 #putting into micro mols
        PAR[dd] = par
        
    return PAR

def Zeucalc(gdata,par):
    zeu = {}
    for dd in ['east','west']:
        PARd = par[dd]
        depth = gdata.zbin.data

        z = np.ones((1,297))*np.nan
        for ix in range(0,297):
            pp = PARd[:,ix]
            nans = np.isnan(pp)
            parTnonan = pp[~nans]
            depthnonan = depth[~nans]
            if len(parTnonan) == 0:
                z[0,ix] = np.nan
            else:
                surf = parTnonan[0]
                try:
                    lims = np.where(pp<surf*0.01)[0]
                    lim = lims[0]
                except IndexError:
                    z[0,ix] = np.nan
                else:
                    z[0,ix] = depth[lim]
        zeu[dd] = z
        
    return zeu

def stratStrengthcalc(gdata):
    stratStrength = {}
    for dd in ['east','west']:
        sigmat = gdata['sigmat_'+dd].data
        surf = np.nanmedian(sigmat[5:10,:],axis=0)
        depth = np.nanmedian(sigmat[57:62,:],axis=0)
        stratStrength[dd] = depth - surf
        
    return stratStrength

def stratStrengthcalc_smooth(gdata):
    #as above, but for smoothed density profiles
    density_smooth = smoothing(gdata,'sigmat',7,3)
    stratStrength = {}
    for dd in ['east','west']:
        sigmat = density_smooth[dd]
        surf = np.nanmedian(sigmat[5:10,:],axis=0)
        depth = np.nanmedian(sigmat[57:62,:],axis=0)
        stratStrength[dd] = depth - surf
        
    return stratStrength

def intChlcalc(gdata):
    intChl = {}
    for dd in ['east','west']:
        chl = gdata['chlfl_'+dd].data
        intChl[dd] = np.nansum(chl[0:60,:],axis=0)
        
    return intChl

def intChlcalc_corr(gdata):
    #as above, but using the quenching
    #corrected chl data
    intChl = {}
    for dd in ['east','west']:
        chl = gdata['chlflCorr_'+dd].data
        intChl[dd] = np.nansum(chl[0:60,:],axis=0)
        
    return intChl

def intChlcalc_corr_depth(gdata,depth):
    #as above, but using the quenching
    #corrected chl data and for input depth
    intChl = {}
    for dd in ['east','west']:
        chl = gdata['chlflCorr_'+dd].data
        intChl[dd] = np.nansum(chl[0:depth,:],axis=0)
        
    return intChl

def avChlcalc(gdata):
    avChl = {}
    for dd in ['east','west']:
        chl = gdata['chlfl_'+dd].data
        avChl[dd] = np.nanmedian(chl[0:65,:],axis=0)
        
    return avChl

def avChlcalc_corr(gdata):
    #as above, but using the quenching
    #corrected chl data
    avChl = {}
    for dd in ['east','west']:
        chl = gdata['chlflCorr_'+dd].data
        avChl[dd] = np.nanmedian(chl[0:65,:],axis=0)
        
    return avChl

def quenchingID(peakData):
    #OLD: was updated after the first paper submission was rejected
    zpeak_idx = (peakData.Zchl < peakData.MLD) & (peakData.Zchl > 0)
    bbC_idx = (peakData.bbC < 20.) & (peakData.bbC > -20.)
    change_idx = (peakData.chlC - peakData.bbC)> 30
    combined_idx = zpeak_idx & bbC_idx & change_idx
    return combined_idx

def QICalc(densityDict,mldDict,depth):
    #quality index for MLD calcs based on Carvalho et al
    #see notebook 01-RevisitingMLD-v2-henry12 for more details
    QI = {}
    directions = ['east','west']
    for dd in directions:
        ddata = densityDict[dd]
        mlddata = mldDict[dd]

        QIdir = np.ones(mlddata.shape)*np.nan
        for ix in range(0,297):
            den = ddata[:,ix]
            try:
                zmld = int(mlddata[ix])
            except ValueError:
                QIdir[ix] = np.nan
            else:
                den1 = den[:zmld]
                den2 = den[:(int(np.round(1.5*zmld)))]
                den1_mean = np.nanmean(den1)
                den2_mean = np.nanmean(den2)

                rmse1 = np.nanstd(den1 - den1_mean,ddof=1)
                rmse2 = np.nanstd(den2 - den2_mean,ddof=1)

                QIdir[ix] = 1 - rmse1/rmse2

            QI[dd] = QIdir.squeeze()
    return QI

def MLDCalc(gdata):
    #see the series of notebooks 01-RevistingMLD for details on
    #creation of this method/function
    
    density_smooth = smoothing(gdata,'sigmat',7,3)

    MLD = {}
    direction = ['east','west']
    for dd in direction:
        N2 = (9.8/density_smooth[dd])*dfdx(density_smooth[dd],gdata.zbin.data)

        mld_smooth = np.ones((1,297))*np.nan
        for ix in range(0,297):
            try:
                N2maxInd = np.nanargmax(N2[:,ix])
                N2max = np.nanmax(N2[:,ix])
                if (N2max > 0.0003) & (gdata.zbin.data[N2maxInd] > 5):
                    mld_smooth[0,ix] = gdata.zbin.data[N2maxInd]
                elif (N2max < 0.0003) & (gdata.zbin.data[N2maxInd] > 5):
                    dataind = np.isnan(N2[:,ix])
                    zind = gdata.zbin.data[~dataind]
                    mld_smooth[0,ix] = np.max(zind)
                else:
                    mld_smooth[0,ix] = np.nan
            except ValueError:
                pass
            
            MLD[dd] = mld_smooth.squeeze()
    
    return MLD