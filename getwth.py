#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Load weather data  """

from __future__ import print_function
import argparse
import numpy as np
import pandas as pd

__autor__ = 'Kyosuke Yamamoto (kyon)'
__date__ = '19 Sep 2017'


def load_config(csvpath):
    ''' Load configs in csvfile '''
    config = {}
    f = open(csvpath)
    lines = f.readlines()
    f.close()
    for line in lines:
        if not line[0] == '#':
            break
        else:
            if line[:10] == '#config - ':
                _config = line.replace('#config - ', '').replace('\n', '')
                key, value = _config.split(':')
                try:
                    config[key] = float(value)
                except ValueError:
                    config[key] = str(value)

    return config


def SVP(temp):
    ''' Saturation vapor pressure (SVP) '''
    return .611 * 10**(7.5 * temp / (237.7 + temp))  # kpa
    #	6.112 * exp(17.67*temp/(243.5 + temp))


def main(csvpath):
    ''' main '''

    #load original
    wth = load_config(csvpath)
    w = pd.read_csv(csvpath, comment='#', sep=',', na_values=['-'])
    w.rename(columns={'YEAR': 'year', 'DOY': 'doy',
                      'swv_dwn': 'srad',
                      'T2M': 'tavg', 'T2MX': 'tmax', 'T2MN': 'tmin',
                      'RAIN': 'prec',
                      'RH2M': 'relh'}, inplace=True)

    #date
    w['date'] = pd.to_datetime(w.year * 1000 + w.doy, format='%Y%j')

    #rhmin, rhmax
    tmin = np.maximum(w.tmin, -5)
    tmax = np.maximum(w.tmax, -5)
    tavg = np.maximum(w.tavg, -5)
    es = SVP(tavg)
    vp = w.relh / 100 * es
    es = SVP(tmax)
    rhmin = 100 * vp / es
    rhmin = np.maximum(0, np.minimum(100, rhmin))
    es = SVP(tmin)
    rhmax = 100 * vp / es
    rhmax = np.maximum(0, np.minimum(100, rhmax))
    w['rhmin'] = rhmin
    w['rhmax'] = rhmax

    #vapr (100 for % and 10 to go from hPa to kPa)
    vapr = w.relh * SVP(w.tavg) / 1000
    w['vapr'] = vapr

    #finalize
    wth['w'] = w

    return wth


if __name__ == '__main__':

    #argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('csvpath', nargs='?', type=str, default='./dataset/daily_weather_28368.nasa.csv')
    args = parser.parse_args()

    wth = main(args.csvpath)
