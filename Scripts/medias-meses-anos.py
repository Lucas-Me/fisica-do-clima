#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 15:07:06 2023

@author: bmiranda
"""
import pandas as pd 

f_latente = pd.read_csv('/home/bmiranda/Desktop/f-clima/trabalho/f-latente-mensal.csv',sep=',')
f_sensivel = pd.read_csv('/home/bmiranda/Desktop/f-clima/trabalho/f-sensivel-mensal.csv',sep=',')

media_meses_latente = f_latente.groupby('mes').mean()
media_meses_sensivel = f_sensivel.groupby('mes').mean()

media_meses_latente.to_csv('media_meses_latente.csv')
media_meses_sensivel.to_csv('media_meses_sensivel.csv')

media_anos_latente = f_latente.groupby('ano').mean()
media_anos_sensivel = f_sensivel.groupby('ano').mean()

media_anos_latente.to_csv('media_anos_latente.csv')
media_anos_sensivel.to_csv('media_anos_sensivel.csv')

