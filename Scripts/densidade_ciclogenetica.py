"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Possui como objetivos criar um arquivo NetCDresults ciclogenetica, a partir
dos ciclones gerados pelo tracking (após o filtro) e um mapa exemplificando os resultados.

Criado por: Lucas Menezes
"""

# IMPORT BUILT-IN MODULES
import os

# IMPORT MODULES
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# IMPORT MAP RELATED MODULES
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def matriz_densidade(
	binx, biny, # Pontos da grade nos eixos
	array_lons : np.array,
	array_lats : np.array
	) -> np.ndarray:
	'''
	Funcao responsavel por calcular a densidade ciclogenetica, mais especificamente
	a quantidade de ciclogeneses em cada pixel da grade.
		
	Argumentos:
		dx = incremento de grade. Ex: 0.5 | Possui unidade de graus 
		array_lons = array com a longitude de cada ponto. (EPSG:4326)
		array_lats = array com a latitude de cada ponto. (EPSG:4326)

		Basicamente ele contabiliza, através dos array de longitude e latitude,
		todos os pontos de um dado elemento, como foco de calor ou raios, que se
		encontram dentro de um pixel (ponto de grade) com resolucao dada por "dx".
	'''

	# Bloco de codigo abaixo faz a contabilizacao da densidade
	H, xedges, yedges = np.histogram2d(array_lons, array_lats, bins = (binx, biny))
	density = H.T  # dimensao 0 (linhas) = latitude, dimensao 1 (colunas)= longitude

	return density


def gerar_netcdf(data_folder, years, resolucao : int = 5):
	'''
	Gera o netcdf com os resultados da densidade ciclogenetica por mes e ano.
	Possui dimensoes (lat, lon, time) e apenas uma variavel (ciclogenese).

	a lista "years" tem que estar em ordem crescente -> [2000, 2001, 2002, ...]
	'''
	
	# cria uma pasta "Processado" no diretorio pai, tudo bem se ja houver
	par_dir = os.path.abspath(os.path.join(data_folder, os.pardir))
	save_folder = os.path.join(par_dir, 'Processado')
	os.makedirs(save_folder, exist_ok=True)

	# resolucao espacial
	dx = resolucao # padrao 5 graus

	# extensao do dominio
	lon_min = -90
	lon_max = -20
	lat_min = -70
	lat_max = -10

	# diumesnoes do netcdf final
	binx = np.arange(lon_min, lon_max + dx, dx) # coordenadas lon
	biny = np.arange(lat_min, lat_max + dx, dx) # coordenadas lat
	times = pd.date_range(start = f"{years[0]}-1-1", end = f"{years[-1]}-12-1", freq = 'MS') # array de datas

	# grade
	x = binx[:-1] + dx / 2
	y = biny[:-1] + dx / 2
	X, Y = np.meshgrid(x, y) # malha da grade
	coords = [y, x, times]
	dims = ['latitude', 'longitude', 'time']
	
	# Valores da variavel
	results = np.zeros((X.shape[0], X.shape[1], times.shape[0])) # começa zerado

	# config
	meses = np.linspace(1, 12, 12).astype(int) # array de meses
	ano_inicio = years[0] # menor ano

	# Loop para cada ano (arquivos):
	for year in years:
		df = pd.read_csv(os.path.join(data_folder, f"{year}.csv"), header = 0)

		# converte lon [0,360] para [-180, 180]
		df['longitude'] = ((df['longitude'] + 180) % 360) - 180

		# convetendo objetos datetime
		df['datetime'] = pd.to_datetime(df['datetime'])

		# loop para cada mes
		for month in meses:
			# slice para o mes e ano
			subset = df.loc[df['datetime'].dt.month == month]

			# resgatando o primeiro passo de tempo de cada ciclone
			subset = subset.loc[subset['step'] == 0]

			# calculando a densidade
			idx = (year - ano_inicio) * 12 + (month - 1)
			results[:, :, idx] =  matriz_densidade(
				binx,
				biny,
				subset['longitude'].to_numpy(),
				subset['latitude'].to_numpy()
				)


	# Gerando o arquivo netcdf
	da = xr.DataArray(results, coords = coords, dims = dims)
	ds = xr.Dataset({"ciclogenese" : da})
	filename = 'frequencia_ciclogenetica.nc'
	ds.to_netcdf(os.path.join(save_folder, filename))


def visualizar_mapas(data_folder, save_folder):
	# abre o arquivo e prepara os dados
	filename = 'frequencia_ciclogenetica.nc'
	ds = xr.open_dataset(os.path.join(data_folder, filename))
	da = ds['ciclogenese']

	# vamo plotar geral.
	total = da.sum(dim = 'time')
	x = total.longitude.values
	y = total.latitude.values
	X, Y = np.meshgrid(x, y)

	# configuracoes do plot
	proj = ccrs.PlateCarree()
	figsize = (10, 10)
	fig, ax  = plt.subplots(
		figsize = figsize,
		subplot_kw= {'projection': proj}
	)

	# Plotando Land e Countries
	ax.add_feature(cf.COASTLINE, color = "lightgray")
	ax.add_feature(cf.BORDERS, color = 'lightgray')

	# Plot points
	CS = ax.contour(
		X, Y,
		total.values
	)
	plt.clabel(CS, inline=1, fontsize=10)
	
	plt.legend(loc = 'upper right')

	# Gridlines
	gl = ax.gridlines(crs=proj, linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
	gl.top_labels = False
	gl.right_labels = False
	gl.xlines = True
	gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 5))
	gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 5))
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

	# extensao
	ax.set_extent([-90, 0, -70, 20], proj)

	# Eixos X e Y
	ax.set_xlabel("LONGITUDE")
	ax.set_ylabel("LATITUDE")

	# mostra o plot
	name = "Total Ciclogeneses.png"
	fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)
