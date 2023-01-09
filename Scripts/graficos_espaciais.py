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
import matplotlib.patches as patches

# IMPORT MAP RELATED MODULES
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def mapa_correlacao(data_folder, save_folder, trimestre = ''):

	suffix = ''
	if len(trimestre) > 0:
		suffix = '_' + trimestre

	# PROPRIEDADES DO PLOT
	# /////////////////////////////////////////////////////
	areas = ['A1', 'A2', 'A3']
	N = len(areas)
	fluxos = {'mslhf' : 'Anomalia FCL', 'msshf' : 'Anomalia FCS'}

	proj = ccrs.PlateCarree()
	fig, axes = plt.subplots(len(fluxos), N, figsize = (12,10), subplot_kw= {'projection' : proj})

	# fronteiras dos paises
	countries = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='50m',
        facecolor='none')

	# LEGENDA 
	# //////////////////////////////////////////////////////

	# define os intervalos da legenda
	levs = np.array(np.linspace(-0.2, .2, 11))

	# lista de cores, em ordem crescete. RGBA
	cmap = mcm.get_cmap('RdBu_r')

	colors = [cmap(i) for i in np.linspace(0, 1, levs.shape[0])]

	# cria um novo cmap
	cmap = mcolors.LinearSegmentedColormap.from_list(
	'Custom cmap', colors, levs.shape[0] - 1)
	#
	cmap.set_under('darkviolet')
	cmap.set_over('darkorange')
	cmap.set_bad('white')

	# nromaliza com base nos intervalos
	norm = mcolors.BoundaryNorm(levs, cmap.N) # usa no PColormesh, nao no Contourf

	# PLOTANDO
	# /////////////////////////////////////////////////////////

	# loop para cada area
	for col in range(N):
		# lendo arquivo
		filename = f'correlacao_mensal_{areas[col]}{suffix}.nc'
		ds = xr.open_dataset(os.path.join(data_folder, filename))

		# coordenadas
		x = ds.longitude.values
		y = ds.latitude.values
		X, Y = np.meshgrid(x, y)

		extent = [x.min(), x.max(), y.min(), y.max()]
		
		# loop para cada fluxo
		row = 0
		for k in fluxos.keys():
			da = ds[k]

			# extensao do mapa
			axes[row, col].set_extent(extent)

			# plota fronteiras
			axes[row, col].add_feature(countries, edgecolor = "k", facecolor= 'none', linewidth = 1)

			# gridlines
			gl = axes[row, col].gridlines(crs=proj, linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
			gl.top_labels = False
			gl.right_labels = False
			gl.xlines = True
			gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 5))
			gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 5))
			gl.xformatter = LONGITUDE_FORMATTER
			gl.yformatter = LATITUDE_FORMATTER

			# Plota os resultados
			axes[row, col].pcolormesh(
				X, Y,
				da.values,
				# levels = levs,
				cmap = cmap,
				norm = norm
				# extend = 'both'
			)

			# titulo
			if row == 0:
				axes[row, col].set_title(f"Região {areas[col]}", fontsize = 20)

			row += 1

	# adiciona legenda 
	fig.subplots_adjust(.1, .1, .9, .9, .2, .1)
	cax = fig.add_axes([.1, .04, .8, .04]) # criando eixo pra colorbar
	cb = fig.colorbar(mcm.ScalarMappable(norm=norm, cmap=cmap), orientation = 'horizontal', pad=0.04, fraction=0.04, cax = cax, extend = 'both')
	font_size = 15 # Adjust as appropriate.
	cb.ax.tick_params(
		labelsize=font_size
		)
	cb.ax.locator_params(nbins=levs.shape[0])

	# labels
	fig.text(0.05, 0.7, "Anomalia FCL", fontsize = 20, ha = 'center', va = 'center', rotation = 'vertical')
	fig.text(0.05, 0.3, "Anomalia FCS", fontsize = 20, ha = 'center', va = 'center', rotation = 'vertical')
	fig.text(0.5, 0.92, f"Coeficiente de Correlação{suffix}", ha = 'center', fontsize = 25)

	# salva a figura e fecha
	filename = f'Correlacao Fluxos e Ciclogeneses{suffix}.png'
	fig.savefig(os.path.join(save_folder, filename), bbox_inches ='tight', dpi = 200)
	plt.close()


def mapa_correlacao_sazonal(data_folder, save_folder):
	trimestres = ['DJF', 'MAM', 'JJA', 'SON']
	for x in range(len(trimestres)):
		mapa_correlacao(data_folder, save_folder, trimestre =  trimestres[x])


def visualizar_mapas(data_folder, save_folder):
	# abre o arquivo e prepara os dados
	filename = 'frequencia_ciclogenetica.nc'
	ds = xr.open_dataset(os.path.join(data_folder, filename))
	da = ds['ciclogenese']

	# vamo plotar geral.
	total = da.sum(dim = 'time')
	total = total * 1e4 # compensacao

	# coordenadas
	x = ds.longitude.values
	y = ds.latitude.values
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
	cmap.set_over('darkred')
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
	ax.set_extent([-80, -30, -60, -15], proj)

	# Eixos X e Y
	ax.set_xlabel("LONGITUDE")
	ax.set_ylabel("LATITUDE")

	ax.set_title("Densidade de ciclogênese [total / km² * 1e-4] (1991 - 2021) ", ha = 'center', fontsize = 15)

	# adicionar areas ciclogenéticas
	rects = [
		patches.Rectangle((-60, -37.5), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj), # A1
		patches.Rectangle((-62.5, -50), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj), # A2
		patches.Rectangle((-72.5, -55), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj) # A3
	]
	leg_rects = ['A1', 'A2', 'A3']
	for i in range(len(rects)):
		x0, y0 = rects[i].get_x(), rects[i].get_y()
		dx, dy = rects[i].get_width(), rects[i].get_height()
		ax.text(x0 + dx / 2, y0 + dy / 2, leg_rects[i], ha = 'center', va = 'center', fontsize = 17, fontweight = 'bold')
		ax.add_patch(rects[i])

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
	q_levs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

	# lista de cores, em ordem crescete. RGBA
	cmap = mcm.get_cmap('Oranges')

	colors = [cmap(i) for i in np.linspace(0, 0.9, q_levs.shape[0])]

	# cria um novo cmap a partir do pre-existente
	cmap = mcolors.LinearSegmentedColormap.from_list(
	'Custom cmap', colors, q_levs.shape[0] - 1)
	# cmap.set_over(np.array([0, 37, 89, 255])/255)
	cmap.set_under('white')
	cmap.set_over('darkred')
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
		axes[lin, col].set_extent([-80, -30, -60, -15], proj)

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
	font_size = 15 # Adjust as appropriate.
	cb.ax.tick_params(labelsize=font_size)

	# salva a figura
	name = "Distribuicao Sazonal Ciclogenese.png"
	fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)


def mapas_tendencia_fluxos(data_folder, save_folder):
	# abre o arquivo e prepara os dados
	filename = 'tendencias_fluxos_1991_2021.nc'
	ds = xr.open_dataset(os.path.join(data_folder, filename))

	vars_ = ['mslhf', 'mslhf anom', 'msshf', 'msshf anom']

	# coordenadas
	x = ds.longitude.values
	y = ds.latitude.values
	X, Y = np.meshgrid(x, y)

	# lendo arquivo dos paises
	countries = cf.NaturalEarthFeature(
		category='cultural',
		name='admin_0_countries',
		scale='50m',
		facecolor='none')
	
	# propriedades dos titulos
	anomalia = {
		'mslhf': 'FCL',
		'mslhf anom' : 'Anomalia Padronizada do FCL',
		'msshf' : 'FCS',
		'msshf anom' : 'Anomalia Padronizada do FCS'
		}

	# loop para cada variavel
	for var in vars_:

		# manipulando
		trend = ds[f'trend {var}'] * 12 # converte de tendencia por mes para tendencia por ano
		trend = trend * 1e2 # compensa os valores mt pequenos

		# mascara onde a condicao é False, ou seja, sem significancia estatistica
		trend = trend.where(ds[f'test {var}']) 

		# LEGENDA
		# /////////////////////////////////////////////////////////////////////////
		
		# percentl para se basear nos intervalos da legenda
		per = np.nanpercentile(trend.values, 1)
		per = np.round(per, 0) # arredondando para a casa mais proxima

		# intervalos da legenda
		q_levs = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

		# lista de cores, em ordem crescete. RGBA
		cmap = mcm.get_cmap('RdBu_r')

		colors = [cmap(i) for i in np.linspace(0, 1, q_levs.shape[0])]

		# cria um novo cmap
		cmap = mcolors.LinearSegmentedColormap.from_list(
		'Custom cmap', colors, q_levs.shape[0] - 1)
		#
		cmap.set_under('darkviolet')
		cmap.set_over('darkorange')
		cmap.set_bad('white')

		# nromaliza com base nos intervalos
		norm = mcolors.BoundaryNorm(q_levs, cmap.N) # usa no PColormesh, nao no Contourf

		# ///////////////////////////////////////////////////////////////////
		# configuracoes do plot
		proj = ccrs.PlateCarree()
		figsize = (10, 10)
		fig, ax  = plt.subplots(
			figsize = figsize,
			subplot_kw= {'projection': proj}
		)

		# Plotando fronteiras dos paises
		ax.add_feature(countries, edgecolor = "gray", facecolor= 'none', linewidth = 1)

		# Plot Contoruf
		CS = ax.contourf(
			X, Y,
			trend.values,
			levels = q_levs,
			cmap = cmap,
			norm = norm,
			extend = 'both'
		)

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
		ax.set_extent([-80, -30, -55, -15], proj)

		# Eixos X e Y
		ax.set_xlabel("LONGITUDE")
		ax.set_ylabel("LATITUDE")

		ax.set_title(f"Tendência da {anomalia[var]} [1/ano x 1e-2]", ha = 'center', fontsize = 15)

		# adicionar areas ciclogenéticas
		rects = [
			patches.Rectangle((-60, -37.5), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj), # A1
			patches.Rectangle((-62.5, -50), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj), # A2
			patches.Rectangle((-72.5, -55), 10, 10, linewidth=3, edgecolor='k', facecolor='none', transform = proj) # A3
		]
		leg_rects = ['A1', 'A2', 'A3']
		for i in range(len(rects)):
			ax.add_patch(rects[i])

		# mostra o plot
		name = f"Tendencia {var}.png"
		fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)
