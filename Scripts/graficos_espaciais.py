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
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.cm as mcm

# IMPORT MAP RELATED MODULES
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def visualizar_mapas(data_folder, save_folder):
	# abre o arquivo e prepara os dados
	filename = 'frequencia_ciclogenetica.nc'
	ds = xr.open_dataset(os.path.join(data_folder, filename))
	da = ds['ciclogenese']

	# vamo plotar geral.
	total = da.sum(dim = 'time')
	total = total * 1e4 # compensacao

	# coordenadas
	x = total.longitude.values
	y = total.latitude.values
	X, Y = np.meshgrid(x, y)
	
	# LEGENDA
	# /////////////////////////////////////////////////////////////////////////
	# define os intervalos da legenda
	q_levs = np.array([1, 2, 3, 4, 5, 7, 10, 13, 17, 25])

	# lista de cores, em ordem crescete. RGBA
	cmap = mcm.get_cmap('Oranges')

	colors = [cmap(i) for i in np.linspace(0, 0.8, q_levs.shape[0])]

	# cria um novo cmap
	cmap = mcolors.LinearSegmentedColormap.from_list(
	'Custom cmap', colors, q_levs.shape[0] - 1)
	#
	cmap.set_under('white')
	cmap.set_over('teal')
	cmap.set_bad('white')

	# nromaliza com base nos intervalos
	norm = mcolors.BoundaryNorm(q_levs, cmap.N) # usa no PColormesh, nao no Contourf
	# norm = mcolors.PowerNorm(1, total.max())

	# configuracoes do plot
	proj = ccrs.PlateCarree()
	figsize = (10, 10)
	fig, ax  = plt.subplots(
		figsize = figsize,
		subplot_kw= {'projection': proj}
	)

	# Plotando fronteiras dos paises
	countries = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='50m',
        facecolor='none')

	ax.add_feature(countries, edgecolor = "k", facecolor= 'none', linewidth = 1)
	# ax.add_feature(cf.BORDERS, color = 'lightgray')

	# Plot Contoruf
	CS = ax.contourf(
		X, Y,
		total.values,
		levels = q_levs,
		cmap = cmap,
		norm = norm,
		extend = 'both'
	)

	# # plot contour
	# ax.contour(
	# 	X, Y,
	# 	total.values,
	# 	levels= q_levs,
	# 	colors = 'k'
	# )

	# adiciona legenda 
	cb = fig.colorbar(CS, orientation = 'horizontal', pad=0.04, fraction=0.04)
	font_size = 10 # Adjust as appropriate.
	cb.ax.tick_params(labelsize=font_size)

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
	ax.set_extent([-90, -30, -70, -15], proj)

	# Eixos X e Y
	ax.set_xlabel("LONGITUDE")
	ax.set_ylabel("LATITUDE")

	ax.set_title("Densidade de ciclogênese [total / km² * 1e-4] (1991 - 2021) ", ha = 'center', fontsize = 15)

	# mostra o plot
	name = "Total Ciclogeneses por km2.png"
	fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)


def sazonalidade_ciclogenese(data_folder, save_folder):
	# abre o arquivo e prepara os dados
	filename = 'frequencia_ciclogenetica.nc'
	ds = xr.open_dataset(os.path.join(data_folder, filename))
	da = ds['ciclogenese']

	# calculando os totais sazonais
	sazonal = da.groupby((da.time.dt.month // 3) % 4).sum()
	sazonal = sazonal * 1e4 # compensacao

	# coordenadas
	x = sazonal.longitude.values
	y = sazonal.latitude.values
	X, Y = np.meshgrid(x, y)
	
	# LEGENDA
	# /////////////////////////////////////////////////////////////////////////
	# define os intervalos da legenda
	q_levs = np.array([1, 3, 5, 7, 9, 11, 13, 15])

	# lista de cores, em ordem crescete. RGBA
	cmap = mcm.get_cmap('Oranges')

	colors = [cmap(i) for i in np.linspace(0, 0.9, q_levs.shape[0])]

	# cria um novo cmap a partir do pre-existente
	cmap = mcolors.LinearSegmentedColormap.from_list(
	'Custom cmap', colors, q_levs.shape[0] - 1)
	# cmap.set_over(np.array([0, 37, 89, 255])/255)
	cmap.set_under('white')
	cmap.set_over(cmap(0.95))
	cmap.set_bad('white')

	# nromaliza com base nos intervalos
	norm = mcolors.BoundaryNorm(q_levs, cmap.N) # usa no PColormesh, nao no Contourf

	# PLOTAGEM
	# ////////////////////////////////////////////////////////////////////////
	
	# configuracoes do plot
	proj = ccrs.PlateCarree()
	figsize = (10, 10)
	trimestre = ['(a) DJF', "(b) MAM", "(c) JJA", "(d) SON"]

	fig, axes  = plt.subplots(
		figsize = figsize,
		subplot_kw= {'projection': proj},
		nrows = 2,
		ncols = 2,
		sharex=True,
		sharey=True
	)

	# fronteiras dos paises
	countries = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='50m',
        facecolor='none')

	# LOOP PARA CADA ESTACAO DO ANO
	for i in range(len(trimestre)):
		lin = i // 2
		col = i % 2

		axes[lin, col].add_feature(countries, edgecolor = "k", facecolor= 'none', linewidth = 1)

		# Plot Contoruf
		CS = axes[lin, col].contourf(
			X, Y,
			sazonal.sel(dict(month = i)),
			levels = q_levs,
			cmap = cmap,
			norm = norm,
			extend = 'both'
		)

		# # plot contour
		# axes[i].contour(
		# 	X, Y,
		# 	total.values,
		# 	levels= q_levs,
		# 	colors = 'k'
		# )

		# Gridlines
		gl = axes[lin, col].gridlines(crs=proj, linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
		gl.top_labels = False
		gl.right_labels = False
		gl.xlines = True
		gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
		gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 10))
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER
	
		# extensao
		axes[lin, col].set_extent([-90, -30, -70, -15], proj)

		# titulo
		axes[lin, col].set_title(trimestre[i], fontsize = '13', loc = 'left')

	# Subplots adjust
	fig.subplots_adjust(
		bottom = .1, 
		top = .9, 
		left = .1, 
		right = .9,
		hspace = 0.1, 
		wspace = 0.1
		)

	# Eixos X e Y e Titulo
	fig.text(0.4, 0.05, "LONGITUDE", fontsize = 15)
	fig.text(0.05, 0.5, "LATITUDE", fontsize = 15, ha = 'center', va = 'center', rotation = 'vertical')
	fig.text(0.5, 0.95, "Densidade de ciclogênese [total / km² * 1e-4] (1991 - 2021) ", ha = 'center', fontsize = 15)

	# adiciona legenda 
	cax = fig.add_axes([.92, .1, .06, .8]) # criando eixo pra colorbar
	cb = fig.colorbar(CS, orientation = 'vertical', pad=0.04, fraction=0.04, cax = cax)
	font_size = 10 # Adjust as appropriate.
	cb.ax.tick_params(labelsize=font_size)

	# salva a figura
	name = "Distribuicao Sazonal Ciclogenese.png"
	fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)

	