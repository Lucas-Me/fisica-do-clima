"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Gera figuras para auxiliar na tomada de decioes, baseando-se em estatistica.

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


def histograma_densidade(data_folder):
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

	bins = np.linspace(1, 25, 25)
	hist, bins = np.histogram(total.to_numpy(), bins = bins)

	norm = mcolors.LogNorm(1, total.max())
	print([norm(i) for i in range(1, 25)])

	fig, ax = plt.subplots()

	ax.bar(bins[:-1], hist, edgecolor = 'none', color = 'royalblue')
	ax.set_yscale('log')
	plt.show()