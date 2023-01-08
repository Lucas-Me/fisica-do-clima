"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Gera as figuras utilizadas nos resultados do seminário.

Criado por: Lucas Menezes
"""

# IMPORT BUILT-IN MODULES
import os

# IMPORT MODULES
import pandas as pd
import numpy as np
import xarray as xr
import pymannkendall as mk

# IMPORT PLOT RELATED MODULES
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.patches as patches

# IMPORT MAP RELATED MODULES
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def tendencia_fluxos(data_folder):
	'''
	CALCULA A TENDENCIA DOS FLUXOS EM CADA PONTO DE GRADE, UTILIZANDO COMO INPUT
	OS DADOS MENSAIS DE CADA ANO.
	Cada calculo leva cerca de 10 a 15 minutos, portanto o calculo da tendencia de todas
	as variaveis deve levar cerca de 1 hora ou mais.
	'''
	# arquivo e vriaveis
	filename = 'SAO-91to21_f-clima.nc'
	vars_ = ['mslhf', 'msshf']

	# lendo arquivo e filtrando por data
	ds = xr.open_dataset(os.path.join(data_folder, filename))
	ds = ds.sel(dict(expver = 1, time = slice('1991-01-01', '2021-12-01')))

	# dicionario com o dataarrays a serem salvos
	new_ds = {}

	# loop para cada variavel
	for var in vars_:

		da = ds[var]

		# agrupando e calculando anomalia
		gp = da.groupby(da.time.dt.month)
		clim = gp.mean(dim = 'time')
		anomalia = gp - clim

		# agrupando e calculando anomalia padronizada
		gp_anom = anomalia.groupby(anomalia.time.dt.month)
		anom_padronizada = gp_anom / clim

		# tendencias para a anomalia padronizada
		trend_anom = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).slope, anom_padronizada, input_core_dims = [['time']], vectorize = True)
		teste_anom = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).h, anom_padronizada, input_core_dims = [['time']], vectorize = True)

		# tendencias para a variavel pura
		trend = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).slope, da, input_core_dims = [['time']], vectorize = True)
		teste = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).h, da, input_core_dims = [['time']], vectorize = True)

		# adiciona no dicionario
		# 'trend' é o próprio valor da tendência
		# 'teste' é o resultado do teste de significancia a um nivel de alpha = 0.05 (True se tem tendencia, False caso contrario)
		new_ds[f'trend {var} anom'] = trend_anom
		new_ds[f'test {var} anom'] = teste_anom
		new_ds[f'trend {var}'] = trend
		new_ds[f'test {var}'] = teste

	# Salva em um netcdf
	filename = 'tendencias_fluxos_1991_2021.nc'
	new_ds = xr.Dataset(new_ds)
	new_ds.to_netcdf(os.path.join(data_folder, filename))