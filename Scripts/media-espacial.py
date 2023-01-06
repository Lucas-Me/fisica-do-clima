#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 10:18:33 2023

@author: bmiranda
"""
#importando bibliotecas
import metpy.calc as mpcalc
import xarray as xr
import pandas as pd

file_1 = xr.open_dataset('/home/bmiranda/Desktop/f-clima/trabalho/data/SAO-91to21_f-clima.nc').metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

level = 1

#RG1
RG1_lon = slice(-49., -36.)
RG1_lat = slice(-24., -33.)

l1 = []
s1 = []

#RG2
RG2_lon = slice(-59., -46.)
RG2_lat = slice(-34., -43.)

l2 = []
s2 = []

#RG3
RG3_lon = slice(-68., -55.)
RG3_lat = slice(-43., -52.)

l3 = []
s3 = []

tempo = []
mes = []
ano = []

lats = file_1.latitude.sel(latitude=RG1_lon).values
lons = file_1.longitude.sel(longitude=RG1_lat).values

for i in range(len(file_1.variables['time'])):
    
    #RG1
    latente1 = file_1.mslhf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG1_lat, 
        longitude=RG1_lon
        ).metpy.unit_array.squeeze()
    
    sensivel1 = file_1.msshf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG1_lat, 
        longitude=RG1_lon
        ).metpy.unit_array.squeeze()
    
    #RG2
    latente2 = file_1.mslhf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG2_lat, 
        longitude=RG2_lon
        ).metpy.unit_array.squeeze()
    
    sensivel2 = file_1.msshf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG2_lat, 
        longitude=RG2_lon
        ).metpy.unit_array.squeeze()
    
    #RG3
    latente3 = file_1.mslhf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG3_lat, 
        longitude=RG3_lon
        ).metpy.unit_array.squeeze()
    
    sensivel3 = file_1.msshf.metpy.sel(
        time = file_1.time[i],
        expver=level,
        latitude=RG3_lat, 
        longitude=RG3_lon
        ).metpy.unit_array.squeeze()
    
    
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    #formatação do time
    t = ('{}'.format(vtime)).split()
    tempo.insert(len(tempo),t[0])
    
    for x in range(len(tempo)):
        m = tempo[x].split('-')
    mes.insert(len(mes),m[1])
    ano.insert(len(ano),m[0])
         
    # calculando as médias de f calor latente para cada tempo e incluindo em listas
    latente1_mean = latente1.mean()
    l1.insert(len(l1),latente1_mean)
    latente2_mean = latente2.mean()
    l2.insert(len(l2),latente2_mean)
    latente3_mean = latente3.mean()
    l3.insert(len(l3),latente3_mean)

    
    #criando dataframe e salvando em .csv
    df = pd.DataFrame()
    df['time'] = tempo
    df['ano'] = ano
    df['mes'] = mes
    df['l1'] = l1
    df['l2'] = l2
    df['l3'] = l3
    df.to_csv('f-latente-mensal.csv')
    
    #calculando as medias de f calor sensivel
    sensivel1_mean = sensivel1.mean()
    s1.insert(len(s1),sensivel1_mean)
    sensivel2_mean = sensivel2.mean()
    s2.insert(len(s2),sensivel2_mean)
    sensivel3_mean = sensivel3.mean()
    s3.insert(len(s3),sensivel3_mean)
    
    #criando dataframe e salvando em .csv
    df = pd.DataFrame()
    df['time'] = tempo
    df['ano'] = ano
    df['mes'] = mes
    df['s1'] = l1
    df['s2'] = l2
    df['s3'] = l3
    df.to_csv('f-sensivel-mensal.csv')
    
    