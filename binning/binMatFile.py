import numpy as np
import pandas as pd
import xarray as xr
import re
import datetime as dt

def binMatFile(mat,glider):
    # finding the variables as named by Bruce
    variablefile = 'V:/Catherine/Glider/DataFiles/Glider Yo Columns.xlsx'
    variables = pd.read_excel(variablefile,sheet_name='new',header=None,usecols=[4],squeeze=True)
    variables = list(variables.drop(variables.index[0]))

    #finding the required variables based on their default names
    if re.search('h',glider):
        searchKeys = ['sci_ocr*','sci_bbfl2s*','m_present_time','x_corrected_lat','x_corrected_lon','x_measured_depth',
                  'x_sci_sigmat','m_pitch','m_roll','sci_water_temp','x_sci_salinity']
    else:
        if 'sci_suna_nitrate_mg' in mat.keys():
            searchKeys = ['sci_ocr*','sci_bbfl2s*','m_present_time','x_corrected_lat','x_corrected_lon','x_measured_depth',
                  'x_sci_sigmat','m_pitch','m_roll','sci_water_temp','x_sci_salinity','sci_oxy4*',
                          'sci_suna_nitrate_mg']
        else:
            searchKeys = ['sci_ocr*','sci_bbfl2s*','m_present_time','x_corrected_lat','x_corrected_lon','x_measured_depth',
                  'x_sci_sigmat','m_pitch','m_roll','sci_water_temp','x_sci_salinity','sci_oxy4*']

    matchedVariables = []
    for key in searchKeys:
        for matkey in mat.keys():
            M = re.match(key,matkey)
            if M:
                matchedVariables += [matkey]

    #dealing with east and west bound data separately
    missionList = []
    directions = ['east','west']
    for direction in directions:

        dirBound = direction + 'bound'

        #selecting required variables from the <direction>bound matrices
        Data = {}
        for vv in matchedVariables:
            ix = int(mat[vv][0][0])
            Data[vv] = mat[dirBound][:,ix-1]

        Data = pd.DataFrame(Data)

        #doing some QC on the data
        subData = Data.loc[np.isfinite(Data['x_corrected_lat'])].copy()
        idx = ~((subData['sci_water_temp'] == 0) & (subData['x_sci_salinity'] < 0) & (subData['x_sci_salinity'] > -1))
        subData = subData.loc[idx].copy()
        subData.dropna(thresh=19,inplace = True)

        #reformatting the datetime string from unix time to year and decimal day of year
        year = []
        jdays = []
        for unixTime in subData['m_present_time']:
            pythonDT = dt.datetime.utcfromtimestamp(unixTime)
            year += [pythonDT.year]
            jday = pythonDT.timetuple().tm_yday
            fractionalDay = (pythonDT.hour/24.) + (pythonDT.minute/60/24) + (pythonDT.second/60/60/24)
            jdays += [jday+fractionalDay]
        subData['year_'+direction] = year
        subData['day_'+direction] = jdays
        subData.drop(['m_present_time'],axis = 1, inplace = True)

        #specifying the depth and longitude bins
        zbins = np.linspace(0,200,num=201)
        lonbins = np.linspace(-69.78,-66.8,num=int((-66.8+69.78)/0.01))

        zlabels = [str(z) for z in range(len(zbins)-1)]
        lonlabels = [str(i) for i in range(len(lonbins)-1)]

        #doing the binning
        subData['zbin'] = pd.cut(subData['x_measured_depth'].values,zbins,labels = zbins[1:])
        subData['lonbin'] = pd.cut(subData['x_corrected_lon'].values,lonbins,labels = lonbins[1:])
        binned = subData.groupby(['zbin','lonbin']).mean()

        #reformatting into xarray
        binned.drop(['x_measured_depth','x_corrected_lon'],axis = 1, inplace=True)
        if 'sci_oxy4_temp' in mat.keys():
            binned.drop(['sci_oxy4_temp'],axis=1,inplace=True)

        data = xr.Dataset.from_dataframe(binned)

        #renaming the variables
        if re.search('h',glider):
            #print('henry')
            if 'sci_ocr507I_irrad1' in mat.keys():
                varnames = {'sci_ocr507R_rad1' : 'Lu412_'+direction,'sci_ocr507R_rad2' : 'Lu443_'+direction,
                        'sci_ocr507R_rad3' : 'Lu490_'+direction,'sci_ocr507R_rad4' : 'Lu510_'+direction,
                        'sci_ocr507R_rad5' : 'Lu532_'+direction,'sci_ocr507R_rad6' : 'Lu555_'+direction,
                        'sci_ocr507R_rad7' : 'Lu670_'+direction,
                        'sci_ocr507I_irrad1' : 'Ed412_'+direction,'sci_ocr507I_irrad2' : 'Ed443_'+direction,
                        'sci_ocr507I_irrad3' : 'Ed490_'+direction,'sci_ocr507I_irrad4' : 'Ed510_'+direction,
                        'sci_ocr507I_irrad5' : 'Ed532_'+direction,'sci_ocr507I_irrad6' : 'Ed555_'+direction,
                        'sci_ocr507I_irrad7' : 'Ed670_'+direction,'sci_bbfl2s_bb_scaled' : 'bb532_'+direction,
                        'sci_bbfl2s_chlor_scaled' : 'chlfl_'+direction,'sci_bbfl2s_cdom_scaled' : 'cdomfl_'+direction,
                        'x_corrected_lat' : 'lat_'+direction,
                        'sci_water_temp' : 'temp_'+direction, 'x_sci_salinity' : 'sal_'+direction,
                        'x_sci_sigmat' : 'sigmat_'+direction,
                        'm_pitch' : 'pitch_'+direction, 'm_roll' : 'roll_'+direction}
            else:
                varnames = {'sci_ocr507r_rad1' : 'Lu412_'+direction,'sci_ocr507r_rad2' : 'Lu443_'+direction,
                        'sci_ocr507r_rad3' : 'Lu490_'+direction,'sci_ocr507r_rad4' : 'Lu510_'+direction,
                        'sci_ocr507r_rad5' : 'Lu532_'+direction,'sci_ocr507r_rad6' : 'Lu555_'+direction,
                        'sci_ocr507r_rad7' : 'Lu670_'+direction,
                        'sci_ocr507i_irrad1' : 'Ed412_'+direction,'sci_ocr507i_irrad2' : 'Ed443_'+direction,
                        'sci_ocr507i_irrad3' : 'Ed490_'+direction,'sci_ocr507i_irrad4' : 'Ed510_'+direction,
                        'sci_ocr507i_irrad5' : 'Ed532_'+direction,'sci_ocr507i_irrad6' : 'Ed555_'+direction,
                        'sci_ocr507i_irrad7' : 'Ed670_'+direction,'sci_bbfl2s_bb_scaled' : 'bb532_'+direction,
                        'sci_bbfl2s_chlor_scaled' : 'chlfl_'+direction,'sci_bbfl2s_cdom_scaled' : 'cdomfl_'+direction,
                        'x_corrected_lat' : 'lat_'+direction,
                        'sci_water_temp' : 'temp_'+direction, 'x_sci_salinity' : 'sal_'+direction,
                        'x_sci_sigmat' : 'sigmat_'+direction,
                        'm_pitch' : 'pitch_'+direction, 'm_roll' : 'roll_'+direction}
        else:
            #print('grampus')
            if 'sci_suna_nitrate_mg' in mat.keys():
                varnames = {'sci_ocr504r_rad1' : 'L380_'+direction,'sci_ocr504r_rad2' : 'Lu443_'+direction,
                    'sci_ocr504r_rad3' : 'Lu490_'+direction,'sci_ocr504r_rad4' : 'Lu532_'+direction,
                    'sci_ocr504i_irrad1' : 'Ed380_'+direction,'sci_ocr504i_irrad2' : 'Ed443_'+direction,
                    'sci_ocr504i_irrad3' : 'Ed490_'+direction,'sci_ocr504i_irrad4' : 'Ed532_'+direction,
                    'sci_bbfl2s_bb_scaled' : 'bb532_'+direction,
                    'sci_bbfl2s_chlor_scaled' : 'chlfl_'+direction,'sci_bbfl2s_cdom_scaled' : 'cdomfl_'+direction,
                    'x_corrected_lat' : 'lat_'+direction,
                    'sci_water_temp' : 'temp_'+direction, 'x_sci_salinity' : 'sal_'+direction,
                    'x_sci_sigmat' : 'sigmat_'+direction,
                    'm_pitch' : 'pitch_'+direction, 'm_roll' : 'roll_'+direction,
                    'sci_suna_nitrate_mg' : 'nitrate_'+direction,
                    'sci_oxy4_oxygen' : 'o2 conc_'+direction, 'sci_oxy4_saturation' : 'o2 saturation_'+direction}
            else:
                varnames = {'sci_ocr504r_rad1' : 'L380_'+direction,'sci_ocr504r_rad2' : 'Lu443_'+direction,
                    'sci_ocr504r_rad3' : 'Lu490_'+direction,'sci_ocr504r_rad4' : 'Lu532_'+direction,
                    'sci_ocr504i_irrad1' : 'Ed380_'+direction,'sci_ocr504i_irrad2' : 'Ed443_'+direction,
                    'sci_ocr504i_irrad3' : 'Ed490_'+direction,'sci_ocr504i_irrad4' : 'Ed532_'+direction,
                    'sci_bbfl2s_bb_scaled' : 'bb532_'+direction,
                    'sci_bbfl2s_chlor_scaled' : 'chlfl_'+direction,'sci_bbfl2s_cdom_scaled' : 'cdomfl_'+direction,
                    'x_corrected_lat' : 'lat_'+direction,
                    'sci_water_temp' : 'temp_'+direction, 'x_sci_salinity' : 'sal_'+direction,
                    'x_sci_sigmat' : 'sigmat_'+direction,
                    'm_pitch' : 'pitch_'+direction, 'm_roll' : 'roll_'+direction,
                    'sci_oxy4_oxygen' : 'o2 conc_'+direction, 'sci_oxy4_saturation' : 'o2 saturation_'+direction}


        data = data.rename(varnames)

        missionList += [data]

    fullMission = xr.merge(missionList)

    return fullMission
