"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Possui como objetivo realizar o pós-processamento do output (planilhas) gerados pelo algoritmo de tracking.

Criado por: Lucas Menezes
"""

# IMPORT BUILT-IN MODULES
import os

# IMPORT MODULES
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# IMPORT MAP RELATED MODULES
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def figura_tracking(tracking : pd.DataFrame, save_folder : str):
	# configuracoes do plot
	proj = ccrs.PlateCarree()
	figsize = (10, 10)
	fig, ax  = plt.subplots(
		figsize = figsize,
		subplot_kw= {'projection': proj}
	)

	# Plotando Land e Countries
	ax.add_feature(cf.COASTLINE, color = "goldenrod")
	ax.add_feature(cf.BORDERS, color = 'goldenrod')

	# Plot points
	ax.plot(tracking['longitude'], tracking['latitude'], marker = 'o', color = 'tab:blue')
	
	# label para cada ponto
	for i, txt in enumerate(tracking['step']):
		texto = f"t + {txt}" if txt > 0 else "t"
		x = tracking.iloc[i]['longitude']
		y = tracking.iloc[i]['latitude']
		ax.annotate(texto, (x, y))

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
	dates = tracking['datetime'].dt.strftime("%d %b %Y %HZ")
	data = dates.loc[tracking['step'] == 0].iloc[0]
	ax.set_title(f"{tracking['id'].iloc[0]} | {data}")

	# mostra o plot
	name = f"Ciclone N{data}.png"
	fig.savefig(os.path.join(save_folder, name), bbox_inches = 'tight', dpi = 200)


def leitura_arquivos(folder, year):
	data_folder = os.path.join(folder, "Dados",'Tracking', 'Processado')
	save_folder = os.path.join(folder, "Figuras")

	# Lendo o arquivo processado
	df = pd.read_csv(os.path.join(data_folder, f"{year}.csv"), header = 0)

	# convertendo objetos datetime
	df['datetime'] = pd.to_datetime(df.datetime)

	# convertendo longitude [0,360] para [-180, 180]
	df['longitude'] = ((df['longitude'] + 180) % 360) - 180

	# id's unicos
	ids = np.unique(df['id'])
	
	# loop entre cada ciclone
	for id_ in ids:
		cyclone = df.loc[df['id'] == id_].copy()
		figura_tracking(cyclone, save_folder = save_folder)