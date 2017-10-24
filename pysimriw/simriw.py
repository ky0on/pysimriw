#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import os
import hjson
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import getwth

__autor__ = 'Kyosuke Yamamoto (kyon)'
__date__ = '20 Sep 2017'


def load_cultivar_params(cultivar_params_file, cultivar):
    ''' Load cultivar parameters '''

    with open(cultivar_params_file, 'r') as f:
        cultivar_params = hjson.load(f)

    if cultivar not in cultivar_params.keys():
        raise Exception('Unknown cultivar:', cultivar + '.', 'Choose from', cultivar_params.keys())

    return cultivar_params[cultivar]


def daylength(lat, doy):
    # if (class(doy) == "Date" | class(doy) == "character") {
    #     doy <- doyFromDate(doy)
    # }

    _doy = doy.copy()

    if lat > 90 or lat < -90:
        lat = np.nan

    _doy.loc[_doy > 365] = 365  # is this ok? uruudoshi?
    if np.any(_doy < 1) or np.any(_doy > 365):
        raise Exception('_doy must be between 1 and 365')

    P = np.arcsin(0.39795 * np.cos(0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.0086 * (doy - 186)))))
    a = (np.sin(0.8333 * np.pi/180) + np.sin(lat * np.pi/180) * np.sin(P))/(np.cos(lat * np.pi/180) * np.cos(P))
    a[a < -1] = -1
    a[a > 1] = 1
    DL = 24 - (24/np.pi) * np.arccos(a)
    return DL


def main(cultivar, weather, transplant, startday, co2, cultivar_params_file='cultivars.hjson', silent=False):
    ''' '''

    #load cultivar parameters
    cultivar = load_cultivar_params(cultivar_params_file, cultivar)

    #load weather data
    wth = getwth.main(weather)

    #Constants and parameters which may not be cutivar specific
    #Constants related to optcal properties of canopy
    SCAT = 0.25
    RSOL = 0.1
    RCAN = 0.22
    KREF = 0.5
    #Parameters related changes with DVI in radiation conversion efficiency
    CL = 0.001
    TAUC = 0.1
    #Conversion factors from rough to brown rice, and panicle to rough rice
    CVBP = 0.76
    CVPP = 0.9

    #Initial conditions for simulation
    if transplant:
        DVII = 0.15  # transplant
        LAII = 0.05
        DWI = 15
    else:
        DVII = 0     # emergence
        LAII = 0.0025
        DWI = 4

    IFLUG1 = 0
    IFLUG2 = 0
    CDED = 0.
    TCHECK = 0.
    STERP = 0.
    HTSUM = 0.
    HTDAY = 0.
    STHT = 0.0
    ATHHT = cultivar['HIMX']
    ATHLT = cultivar['HIMX']
    JJ = 1
    DWGRAIN = 0.0
    DVI = DVII
    LAI = LAII
    DW = DWI
    DWPAN = 0.0
    STHT = 0
    STLT = 0

    #weather data
    AVT = wth['w'].tavg
    RAD = wth['w'].srad
    TMX = wth['w'].tmax
    startday = pd.to_datetime(startday)
    endday = startday + pd.to_timedelta('200 days')
    days = pd.date_range(startday, endday)
    DL = daylength(wth['lat'], wth['w'].doy)

    startindex = wth['w'].date.eq(startday).idxmax()
    day = startindex - 1
    growing = True
    simday = 0

    res = {}
    # res = as.data.frame(matrix(ncol=9, nrow=length(days)))
    # colnames(res) = c('date','TMP', 'RAD','DL','DVI','LAI', 'DW', 'GY', 'PY')
    # class(res[,'date']) = 'Date'

    #Dynamic Section of The Model
    while growing:
        day += 1
        simday += 1
        if day >= wth['w'].shape[0] - 1:
            warnings.warn('reached end of weather records')
            growing = False

        res[simday] = {
            'date': wth['w'].loc[day, 'date'],
            'TMP': AVT[day],
            'RAD': RAD[day],
            'DL': DL[day],
            'DVI': DVI,
            'LAI': LAI,
            'DW': DW,
            'GY': DWGRAIN,
            'PY': DWPAN,
        }

        #Culculation of Developmental Index DVI
        if DVI < cultivar['DVIA']:
            #before crop becomes photo sensitive period
            EFT = AVT[day] - cultivar['TH']
            DVR = 1. / (cultivar['GV'] * (1.0 + np.exp(-cultivar['ALF'] * EFT)))
        elif DVI <= 1.0:
            #before heading (in photo sensitive period)
            EFT = AVT[day] - cultivar['TH']
            EFL = min(DL[day] - cultivar['LC'], 0.)
            DVR = (1.0 - np.exp(cultivar['BDL'] * EFL)) / (cultivar['GV'] * (1.0 + np.exp(-cultivar['ALF'] * EFT)))
        else:
            #between heading and maturity
            EFT = max(AVT[day] - cultivar['TCR'], 0.)
            DVR = (1.0 - np.exp(-cultivar['KCR'] * EFT)) / cultivar['GR']
        DVI += DVR

        #Culculation of LAI
        if DVI < 0.95:
            #before heading
            EFFTL = max(AVT[day] - cultivar['TCF'], 0.)  # effective temperature for LAI growth (>=0)
            GRLAI = LAI * cultivar['A'] * \
                (1.0 - np.exp(-cultivar['KF'] * EFFTL)) * \
                (1.0 - (LAI/cultivar['FAS'])**cultivar['ETA'])   # growth rate of LAI
            GRL95 = GRLAI
            DVI95 = DVI
        elif GRLAI > 0.0 or DVI <= 1.0:
            #no explanation
            GRLAI = GRL95 * (1.0 - (DVI - DVI95)/(1 - DVI95))
            LAIMX = LAI
            DVI1 = DVI
        elif DVI < 1.1:
            #after heading (empirical function)
            GRLAI = -(LAIMX * (1.0 - cultivar['BETA']) * (DVI - DVI1)/(1.1 - DVI1))*DVR
        else:
            #after heading (empirical function)
            GRLAI = -LAIMX * (1.0 - cultivar['BETA']) * DVR
        LAI += GRLAI

        #Culuculation of Crop Dry Weight
        # SCAT: scattering coefficient
        # RCAN: reflectance when the surface is completely covered by the vegetation
        # RSOL: reflectance of bare soil
        # REF:  canopy reflectance
        # KREF: maybe empirical constant
        # ABSRAD: amount of radiation absorbed by the canopy
        # CONEF, COVCO2: radiation conversion efficiency
        TAU = np.exp(-(1.0 - SCAT) * cultivar['EXTC'] * LAI)
        REF = RCAN - (RCAN - RSOL) * np.exp(-KREF * LAI)   # canopy reflectance, KREF=1/2
        ABSOP = 1.0 - REF - (1.0 - RSOL) * TAU
        ABSRAD = RAD[day] * ABSOP
        COVCO2 = cultivar['COVES'] * (1.54 * (co2 - 330.0)/(1787.0+(co2 - 330.0))+1.0)
        if DVI < 1.0:
            CONEF = COVCO2
        else:
            #decreases gradually toward zero at maturity
            CONEF = COVCO2 * (1.0+CL)/(1.0+CL * np.exp((DVI - 1.0)/TAUC))
        DW = DW + CONEF * ABSRAD  # daily dry matter production = absorbed radiation * radiation conversion efficiency

        #Culuculation of Spikelet Sterility Percentage due to Cool Temerature
        #(the period of highest sensitivity of the rice panicle to cool temperatures)
        # SSTR: Spikelet Sterility bla bla (maybe)
        if DVI > 0.75 and DVI < 1.2:
            CDEG = max(cultivar['THOT'] - AVT[day], 0.)  # effective cool temperature
            CDED += CDEG   # cooling degree-days
            SSTR = cultivar['STO'] + cultivar['BST'] * CDED**cultivar['PST']   # Eq. 4.14
            STLT = min(100.0, SSTR)   # must be <=100
            RIPEP = 1.0 - STLT/100.0  # 0-100 -> 1-0
            ATHLT = cultivar['HIMX'] * RIPEP   # harvest index considering cool temperature

        #Culculation of Spikelet Sterility Percentage due to Heat Stress
        if DVI > 0.96 and DVI < 1.20:
            #accumulate daily max temperature around anthesis
            HTSUM += TMX[day]
            HTDAY += 1

        if DVI >= 1.20 and IFLUG1 == 0:
            AVTMX = HTSUM / HTDAY   # average daily max temperature around anthesis
            STHT = 100.0/(1.0 + np.exp(-0.853 * (AVTMX - 36.6)))  # Eq. 4.16 (Relation between average daily maximum temperature during the floweringperiod and spikelet fertility)
            ATHHT = (1.0 - STHT/100.0) * cultivar['HIMX']   # 0-100 -> 1-0, then harvest index considering high temperature
            IFLUG1 = 1

        #Culculation of Grain Yield
        # HI: harvest index
        ATHI = min(ATHLT, ATHHT)  # choose the largest effect from cool and heat stresses (actual spikelet sterility)
        STERP = max(STHT, STLT)                    # not used
        EFDVI = max(DVI - 1.22, 0.0)               # DVI after flowering
        HI = ATHI * (1.0 - np.exp(-5.57 * EFDVI))  # harvest index
        DWGRAIN = DW * HI
        DWPAN = DWGRAIN/CVBP/CVPP

        #Time Control and Terminal Condition of Simulation
        if DVI > 1.0 and AVT[day] < cultivar['CTR']:
            TCHECK += 1
        if DVI > 2.0:
            if not silent:
                print('DVI reached to 2.')
            growing = False

    #finalize
    simday += 1
    res[simday] = {
        'date': wth['w'].loc[day, 'date'],
        'TMP': AVT[day],
        'RAD': RAD[day],
        'DL': DL[day],
        'DVI': DVI,
        'LAI': LAI,
        'DW': DW,
        'GY': DWGRAIN,
        'PY': DWPAN,
    }

    #Terminal Section of  Simulation
    PYBROD = DWGRAIN / 100.0
    PYBRON = PYBROD / 0.86
    PYPADY = PYBRON / CVBP
    PANDW = PYBROD / CVBP / CVPP
    DWT = DW / 100.0

    #simulation result
    simulated = {}
    simulated['cultivar'] = cultivar
    simulated['PYBROD'] = PYBROD
    simulated['PYBRON'] = PYBRON
    simulated['PYPADY'] = PYPADY
    simulated['PANDW'] = PANDW
    simulated['DWT'] = DWT
    simulated['d'] = pd.DataFrame(res).T   # [1:simday, ]

    return simulated


if __name__ == '__main__':

    #argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--cultivar', '-c', default='Nipponbare', type=str)
    parser.add_argument('--weather', '-w', default='./dataset/daily_weather_28368.nasa.csv', type=str)
    parser.add_argument('--startday', '-s', default='2000-05-15', type=str)
    parser.add_argument('--co2', default=350, type=int)
    parser.add_argument('--transplant', action='store_true')
    parser.add_argument('--out', default='output', type=str)
    args = parser.parse_args()

    #init
    plt.style.use('ggplot')
    simulated = main(args.cultivar, args.weather, args.transplant,
                     args.startday, args.co2,
                     cultivar_params_file='cultivars.hjson')

    #plot
    simulated['d'][['DW', 'GY', 'PY']].plot()
    plt.savefig(os.path.join(args.out, 'simulated.pdf'))
    simulated['d'].to_csv(os.path.join(args.out, 'simulated.csv'))
    print('\nsimulated["d"].tail():\n', simulated['d'].tail())
