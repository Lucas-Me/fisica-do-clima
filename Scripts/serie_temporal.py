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
import pymannkendall as mk

# MODULES PARA PLOT
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm

def ciclogenese_temporal(data_dir, save_folder):

	# lendo o arquivo
	df = pd.read_csv(os.path.join(data_dir, 'frequencia_ciclogenese_mensal_AREAS.csv'))

	# agrupando e aplicando operacoes
	media_mensal = df.groupby(['month']).agg(
		A1 = ('A1', 'mean'),
		A2 = ('A2', 'mean'),
		A3 = ('A3', 'mean')
	)

	total_anual = df.groupby(['ano']).agg(
		A1 = ('A1', 'sum'),
		A2 = ('A2', 'sum'),
		A3 = ('A3', 'sum')
	)

	anom_padronizada = (total_anual - total_anual.mean()) / total_anual.std()
	print(anom_padronizada)

	# figuras
	fig, axes = plt.subplots(nrows = 2, figsize = (15, 10)) # 1 plot mensal e outro anual

	# propriedades
	meses = ['Jan', 'Fev', 'Mar', 'Abr', 'Maio', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dec']
	cmap = mcm.get_cmap('CMRmap') # pega cores a partir de um colormap
	areas = ['A1', 'A2', 'A3']
	cores = cmap(np.linspace(0, 1, len(areas) + 2))

	# Plotando para cada area
	for i in range(len(areas)):
		# plot mensal
		p = axes[0].plot(media_mensal.index, media_mensal[areas[i]], color = cores[i + 1], marker = 'o', markerfacecolor= 'none', label = areas[i], linewidth = 2)

		# plot anual
		axes[1].plot(total_anual.index, total_anual[areas[i]], color = cores[i + 1], marker = 'o', markerfacecolor= 'none', label = areas[i], linewidth = 2)

	# eixo horizontal
	axes[0].set_xticks(media_mensal.index)
	axes[0].set_xticklabels(meses, fontsize = 15)
	axes[0].grid(True, axis = 'y')

	axes[1].set_xticks(np.arange(1991, 2022, 5).astype(int))
	axes[1].set_xticklabels(np.arange(1991, 2022, 5).astype(int), fontsize = 15)
	axes[1].grid(True, axis = 'both')

	# titulos
	axes[0].set_title('(a) Médias anuais', fontsize = 15, loc = 'left')
	axes[1].set_title('(b) Totais anuais', fontsize = 15, loc = 'left')

	# ajustando a figura
	fig.subplots_adjust(0.1, 0.1, 0.9, 0.9, 0.1, 0.2)
	fig.text(0.05, 0.5, 'Frequência de ciclogênese por ano', fontsize = 20, va = 'center', rotation = 'vertical')

	# Legenda
	bbox = [0.5, 0.92]
	handles, labels = axes[0].get_legend_handles_labels()
	lg = fig.legend(handles, labels, ncols = 3, bbox_to_anchor = bbox, loc = 'center', frameon = False, fontsize = 20)

	# mostra a figura
	filename ='serie_temporal_ciclogenese.png'
	fig.savefig(os.path.join(save_folder, filename), bbox_inches = 'tight')
	plt.close()