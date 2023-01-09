"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Gera as figuras utilizadas nos resultados do seminário.

Criado por: Lucas Menezes
"""

# IMPORT BUILT-IN MODULES
import os
from datetime import date

# IMPORT MODULES
import pandas as pd
import numpy as np
import xarray as xr
import pymannkendall as mk
from scipy.stats import spearmanr
from scipy.fft import fft, fftfreq

# IMPORT PLOT RELATED MODULES
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm


def transformada_fourier(data_folder, save_folder):
	# lendo o arquivo
	df = pd.read_csv(os.path.join(data_folder, 'frequencia_ciclogenese_mensal_AREAS.csv'))

	total_anual = df.groupby(['ano']).agg(
		A1 = ('A1', 'sum'),
		A2 = ('A2', 'sum'),
		A3 = ('A3', 'sum')
	)

	# areas ciclogeneticas
	areas = ['A1', 'A2', 'A3']
	
	# plot
	cmap = mcm.get_cmap('CMRmap')
	cores = cmap(np.linspace(0, 1, len(areas) + 2))
	fig, ax = plt.subplots(figsize = (8, 5))

	#loop para cada area
	for i in range(len(areas)):
		N = total_anual.shape[0]
		T = 1
		yf = fft(total_anual[areas[i]].values)
		xf = fftfreq(N, T)[:N//2]
		potencia = 2/N * np.abs(yf[0:N//2])
		xf = xf[1:]
		potencia = potencia[1:]
		anos = 1 / xf
		ax.plot(xf, potencia, color = cores[i + 1], label = areas[i], linewidth = 2)

		# Printando a origem das 3 primeiras frequencias
		ordenado = np.sort(potencia)
		print(ordenado)
		print('\n', areas[i])
		for j in range(3):
			value = xf[np.where(potencia == ordenado[- (j + 1)])[0][0]]
			period = 1 / value
			print(f'{j} || freq = {value:.2f} / ano || periodo {period:.2f} anos')

	# ajustando a figura
	fig.subplots_adjust(0.1, 0.1, 0.9, 0.9, 0.1, 0.1)
	fig.text(0.04, 0.5, 'Potência', fontsize = 15, va = 'center', rotation = 'vertical')
	fig.text(0.5, 0.02, 'Frequência de oscilação [1 / ano]', fontsize = 15, ha = 'center')

	# Legenda
	bbox = [0.5, 0.92]
	handles, labels = ax.get_legend_handles_labels()
	lg = fig.legend(handles, labels, ncols = 3, bbox_to_anchor = bbox, loc = 'center', frameon = False, fontsize = 15)

	ax.grid(True, axis = 'both')
	
	# salva
	filename = 'transformada_fourier_anual_ciclogenese.png'
	fig.savefig(os.path.join(save_folder, filename), bbox_inches = 'tight', dpi = 200)
	plt.close()


def correlacao_fluxos_ciclogenese_sazonal(era5_dir, ciclo_dir):
	'''
	CALCULA A CORRELACAO DOS FLUXOS EM CADA PONTO DE GRADE, UTILIZANDO COMO INPUT
	OS DADOS MENSAIS DE ANOMALIA PADRONIZADA DE CADA ANO.
	'''
	def correlation(x1, x2, alpha):
		
		result = spearmanr(x1, x2)
		corr = result.correlation
		if result.pvalue >= alpha:
			return np.nan

		return corr # to return a single correlation index, instead of a matriz

	# manipulacao dos dados
	# /////////////////////////////////////////////////////////////////////

	# arquivo e vriaveis
	filename = 'SAO-91to21_f-clima.nc'

	# lendo arquivo e filtrando por data
	ds = xr.open_dataset(os.path.join(era5_dir, filename))
	ds = ds.drop_vars(['msl'], errors ='ignore')
	ds = ds.sel(dict(expver = 1, time = slice('1991-01-01', '2021-12-01')))

	def seasonal_years_split(dates):
		years = dates.dt.year
		seasons = ((dates.dt.month // 3) % 4).values
		months = ['01', '04', '07', '10']
		new_dates = [date.fromisoformat(f'{years[i].values}-{months[seasons[i]]}-01') for i in range(seasons.shape[0])]
		new_dates = np.array(new_dates, dtype = np.datetime64)
		new_dates = xr.DataArray(name = 'time', data = new_dates, dims= dates.dims, coords = dates.coords)

		return new_dates

	# media por estacao e ano
	split = seasonal_years_split(ds.time)
	ds = ds.groupby(split).mean(dim = 'time')

	# Calculando as anomalias
	gp = ds.groupby(ds.time.dt.month)
	clim = gp.mean(dim = 'time')
	anomalia = gp - clim

	# agrupando e calculando anomalia padronizada
	gp_anom = anomalia.groupby(anomalia.time.dt.month)
	anom_padronizada = gp_anom / gp.std(dim = 'time')

	# lendo o arquivo da ciclogenese
	df = pd.read_csv(os.path.join(ciclo_dir, 'frequencia_ciclogenese_mensal_AREAS.csv'))
	df['season'] = (df.month // 3) % 4
	df = df.groupby(['season', 'ano']).agg(
		A1 = ('A1', 'sum'),
		A2 = ('A2', 'sum'),
		A3 = ('A3', 'sum'),
	).reset_index(drop = False)

	# slices para cada area, em ordem A1, A2, A3
	areas = ['A1', 'A2', 'A3']
	slice_lon = [slice(-60.0, -50.0), slice(-62.5, -52.5), slice(-72.5, -62.5)]
	slice_lat = [slice(-27.5, -37.5), slice(-40.0, -50.0), slice(-45.0, -55.0)]
	
	# anomalias para a ciclogense
	group = df.groupby(df['season'])
	anom = df[areas] - group[areas].transform('mean')
	std_anom = anom / group[areas].transform('std')
	std_anom['season'] = df['season']
	print(std_anom)

	# Configuracao do plot
	# //////////////////////////////////////////////////////////////////////

	estacao = ['DJF', 'MAM', 'JJA', 'SON']
	for j in range(4):
		season_ds = anom_padronizada.sel(dict(time = (anom_padronizada.time.dt.month // 3) % 4 == j))
		this_df = std_anom.loc[std_anom['season'] == j]
		for i in range(len(areas)):
			area = season_ds.sel(dict(longitude = slice_lon[i], latitude = slice_lat[i]))
			corr_area = xr.apply_ufunc(correlation, area, input_core_dims=[['time']], vectorize=True, kwargs = dict(x2 = this_df[areas[i]], alpha = 0.05))

			# salva 
			filename = f'correlacao_mensal_{areas[i]}_{estacao[j]}.nc'
			corr_area.to_netcdf(os.path.join(era5_dir, filename))
		

def correlacao_fluxos_ciclogenese(era5_dir, ciclo_dir):
	'''
	CALCULA A CORRELACAO DOS FLUXOS EM CADA PONTO DE GRADE, UTILIZANDO COMO INPUT
	OS DADOS MENSAIS DE ANOMALIA PADRONIZADA DE CADA ANO.
	'''
	def correlation(x1, x2, alpha):
		
		result = spearmanr(x1, x2)
		corr = result.correlation
		if result.pvalue >= alpha:
			return np.nan

		return corr # to return a single correlation index, instead of a matriz

	# manipulacao dos dados
	# /////////////////////////////////////////////////////////////////////

	# arquivo e vriaveis
	filename = 'SAO-91to21_f-clima.nc'
	vars_ = ['mslhf', 'msshf']

	# lendo arquivo e filtrando por data
	ds = xr.open_dataset(os.path.join(era5_dir, filename))
	ds = ds.drop_vars(['msl'], errors ='ignore')
	ds = ds.sel(dict(expver = 1, time = slice('1991-01-01', '2021-12-01')))

	# Calculando as anomalias
	gp = ds.groupby(ds.time.dt.month)
	clim = gp.mean(dim = 'time')
	anomalia = gp - clim

	# agrupando e calculando anomalia padronizada
	gp_anom = anomalia.groupby(anomalia.time.dt.month)
	anom_padronizada = gp_anom / gp.std(dim = 'time')

	# lendo o arquivo da ciclogenese
	df = pd.read_csv(os.path.join(ciclo_dir, 'frequencia_ciclogenese_mensal_AREAS.csv'))
	
	areas = ['A1', 'A2', 'A3']
	# slices para cada area, em ordem A1, A2, A3
	slice_lon = [slice(-60.0, -50.0), slice(-62.5, -52.5), slice(-72.5, -62.5)]
	slice_lat = [slice(-27.5, -37.5), slice(-40.0, -50.0), slice(-45.0, -55.0)]
	
	# anomalias para a ciclogense
	group = df.groupby(df['month'])
	anom = df[areas] - group[areas].transform('mean')
	std_anom = anom / group[areas].transform('std')

	# Configuracao do plot
	# //////////////////////////////////////////////////////////////////////

	for i in range(len(areas)):
		area = anom_padronizada.sel(dict(longitude = slice_lon[i], latitude = slice_lat[i]))
		corr_area = xr.apply_ufunc(correlation, area, input_core_dims=[['time']], vectorize=True, kwargs = dict(x2 = std_anom[areas[i]], alpha = 0.05))

		# salva 
		filename = f'correlacao_mensal_{areas[i]}.nc'
		corr_area.to_netcdf(os.path.join(era5_dir, filename))


def tendencia_fluxos(data_folder):
	'''
	CALCULA A TENDENCIA DOS FLUXOS EM CADA PONTO DE GRADE, UTILIZANDO COMO INPUT
	OS DADOS MENSAIS DE CADA ANO.
	Cada calculo leva cerca de 10 a 15 minutos, portanto o calculo da tendencia de todas
	as variaveis deve levar cerca de 2 hora ou mais.
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
		print('Iniciando os procedimentos para a variavel', var)

		da = ds[var]

		# agrupando e calculando anomalia
		gp = da.groupby(da.time.dt.month)
		clim = gp.mean(dim = 'time')
		anomalia = gp - clim

		# agrupando e calculando anomalia padronizada
		gp_anom = anomalia.groupby(anomalia.time.dt.month)
		anom_padronizada = gp_anom / gp.std(dim = 'time')

		# tendencias para a anomalia padronizada
		print(f"Calculando a tendencia de {var} (anomalia)")
		trend_anom = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).slope, anom_padronizada, input_core_dims = [['time']], vectorize = True)

		print(f"Calculando o teste de significancia de {var} (anomalia)")
		teste_anom = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).h, anom_padronizada, input_core_dims = [['time']], vectorize = True)

		# tendencias para a variavel pura
		print(f"Calculando a tendencia de {var}")
		trend = xr.apply_ufunc(lambda x: mk.yue_wang_modification_test(x).slope, da, input_core_dims = [['time']], vectorize = True)

		print(f"Calculando o teste de significancia de {var}")
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