#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import argparse
import matplotlib.pyplot as plt

from pysimriw import simriw

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
    simulated = simriw.main(args.cultivar, args.weather, args.transplant,
                            args.startday, args.co2,
                            cultivar_params_file='cultivars.hjson')

    #plot
    simulated['d'][['DW', 'GY', 'PY']].plot()
    plt.savefig(os.path.join(args.out, 'simulated.pdf'))
    simulated['d'].to_csv(os.path.join(args.out, 'simulated.csv'))
    print('\nsimulated["d"].tail():\n', simulated['d'].tail())
