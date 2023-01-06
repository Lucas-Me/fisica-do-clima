"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de
Física do Clima.

Possui como objetivo criar figuras com os intervalos interquartis da vorticidade
relativa e sua derivada temporal em cada passo de tempo dos ciclones.

O tracking dos ciclones foi avaliado para o ano de 2021 utilizando cartas sinóticas
da marinha, e assim classificados de acordo com a coluna "VALIDO" em:

LEGENDA DA COLUNA "VALIDO"

4 -> CICLONE NOVO QUE SE DEU A PARTIR DO RAMO FRONTAL DE UM EXTRATROPICAL ANTIGO
3 -> CICLONE EM ESTAGIO MADURO / OCLUSAO
2 -> CICLONE NOVO QUE SE DEU A PARTIR DA DIVISAO DE UM EXTRATROPICAL ANTIGO
1 -> CICLOGENESE
0 -> ERRADO (nao se trata de um ciclone, PODE ser um cavado)

A intenção é utilizar os dados validados para encontrar um padrão nos dados de ciclogenese
e não ciclogenese, a fim de distinguir entre estes dois. Este padrao sera utilizado
mais tarde para calibração do tracking nos demais anos (aplicação de um filtro).

Criado por: Lucas Menezes
"""

# IMPORT MODULES
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# PROPRIEDADES DO MATPLOTLIB
rcparams = {
	'font.family' : 'Times New Roman'
}
plt.rcParams.update(rcparams)

def interquartil_vorticidade_derivada(df : pd.DataFrame, save_folder):

	# Intervalos interquartis nas primeiras 48h
	group = df.groupby(['legenda', 'step']) # agrupando por passo de tempo

	# coluna legenda
	data_legend = np.unique(df['legenda'])
	
	# resultados de percentil e média
	results = group.agg(
		P5 = ('derivada_1', lambda x: np.nanpercentile(x, 5)), # percentil 5
		P25 = ('derivada_1', lambda x: np.nanpercentile(x, 25)), # percenntil 25
		P50 = ('derivada_1', lambda x: np.nanpercentile(x, 50)), # mediana
		P75 = ('derivada_1', lambda x: np.nanpercentile(x, 75)), # percentil 75
		P95 = ('derivada_1', lambda x: np.nanpercentile(x, 95)), # percentil 95
		MEDIA = ('derivada_1', 'mean') # media
	)

	# Configuracoes do grafico
	fig, axes = plt.subplots(figsize = (12, 12), nrows = data_legend.shape[0], sharex = True)

	# Plotando os resultados
	for i in range(data_legend.shape[0]):
		subset = results.loc[(data_legend[i],)]

		# Mediana  e media
		axes[i].plot(subset.index, subset['P50'], color = 'green', linewidth = 2)
		# axes[i].plot(subset.index, subset['MEDIA'], color = 'royalblue', linewidth = 2)

		# preenche entre o percentil 25 e 75
		axes[i].fill_between(subset.index, subset['P25'], subset['P75'], color = 'orangered')

		# preenche entre o percentil 5 e 25
		axes[i].fill_between(subset.index, subset['P5'], subset['P25'], color = 'orange')

		# preenche entre o percentil 75 e 95
		axes[i].fill_between(subset.index, subset['P75'], subset['P95'], color = 'orange')

		# FFORMATACAO DO EIXO X.
		axes[i].set_xticks(subset.index)
		axes[i].set_xticklabels(["t"] + [f"t + {passo * 6}h" for passo in subset.index[1:]], fontsize = 10)
		axes[i].set_xlim(0, 14) # Primeiras 48 horas

		# FORMATACOES DO EIXO Y
		axes[i].set_ylim(-0.3, 0.5)
		axes[i].set_yticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5])

		# Titulo
		axes[i].set_title(f"CICLONES {data_legend[i]} (2021)", fontsize = 13)

		# Demais configuracoes
		axes[i].grid(True, axis = 'y')


	# Legenda do eixo x
	fig.text(0.5, 0.07, "Passo de tempo", ha = 'center', fontsize = 18)

	# Legenda do eixo y
	fig.text(0.07, 0.5, "Tendencia da vorticidade relativa", va = 'center', rotation = 'vertical', fontsize = 18)

	# Salva a figura e fecha
	nome = os.path.join(save_folder, "Intervalos Interquartis Derivada.png")
	fig.savefig(nome, bbox_inches = 'tight', dpi = 200)
	plt.close()


def interquartil_vorticidade(df : pd.DataFrame, save_folder):

	# Intervalos interquartis nas primeiras 48h
	group = df.groupby(['legenda', 'step']) # agrupando por passo de tempo

	# coluna legenda
	data_legend = np.unique(df['legenda'])
	
	# resultados de percentil e média
	results = group.agg(
		P5 = ('rel_vort', lambda x: np.nanpercentile(x, 5)), # percentil 5
		P25 = ('rel_vort', lambda x: np.nanpercentile(x, 25)), # percenntil 25
		P50 = ('rel_vort', lambda x: np.nanpercentile(x, 50)), # mediana
		P75 = ('rel_vort', lambda x: np.nanpercentile(x, 75)), # percentil 75
		P95 = ('rel_vort', lambda x: np.nanpercentile(x, 95)), # percentil 95
		MEDIA = ('rel_vort', 'mean') # media
	)

	# Configuracoes do grafico
	fig, axes = plt.subplots(figsize = (12, 12), nrows = data_legend.shape[0], sharex = True)

	# Plotando os resultados
	for i in range(data_legend.shape[0]):
		subset = results.loc[(data_legend[i],)]

		# Mediana  e media
		axes[i].plot(subset.index, subset['P50'], color = 'green', linewidth = 2)
		# axes[i].plot(subset.index, subset['MEDIA'], color = 'royalblue', linewidth = 2)

		# preenche entre o percentil 25 e 75
		axes[i].fill_between(subset.index, subset['P25'], subset['P75'], color = 'orangered')

		# preenche entre o percentil 5 e 25
		axes[i].fill_between(subset.index, subset['P5'], subset['P25'], color = 'orange')

		# preenche entre o percentil 75 e 95
		axes[i].fill_between(subset.index, subset['P75'], subset['P95'], color = 'orange')

		# FFORMATACAO DO EIXO X.
		axes[i].set_xticks(subset.index)
		axes[i].set_xticklabels(["t"] + [f"t + {passo * 6}h" for passo in subset.index[1:]], fontsize = 10)
		axes[i].set_xlim(0, 14) # Primeiras 48 horas

		# FORMATACOES DO EIXO Y
		axes[i].set_ylim(0, 10)

		# Titulo
		axes[i].set_title(f"CICLONES {data_legend[i]} (2021)", fontsize = 13)

		# Demais configuracoes
		axes[i].grid(True, axis = 'y')


	# Legenda do eixo x
	fig.text(0.5, 0.07, "Passo de tempo", ha = 'center', fontsize = 18)

	# Legenda do eixo y
	fig.text(0.07, 0.5, "Vorticidade relativa negativa", va = 'center', rotation = 'vertical', fontsize = 18)

	# Salva a figura e fecha
	nome = os.path.join(save_folder, "Intervalos Interquartis.png")
	fig.savefig(nome, bbox_inches = 'tight', dpi = 200)
	plt.close()


def sitemas_semiestacionarios(tracking : pd.DataFrame):
	'''
	Identifica sistemas semiestacionarios de acordo com o seu deslocamento total.
	O criterio utilizado para essa identificacao foi possuir deslocamento inferior
	a 10°.
	'''
	ids = np.unique(tracking['id'])
	semi_estacionario = pd.Series(False, tracking.index)

	for id_ in ids:
		condition = tracking['id'] == id_
		cyclone = tracking.loc[condition]
		
		# ordena de acordo com o passo de tempo
		cyclone = cyclone.sort_values(by = 'step', ascending=True)

		distance = 0 # distancia total em graus
		for row in cyclone.iterrows():
			row = row[1]
			if row['step'] == 0: # se for o primeiro passo de tempo, pula.
				continue
			
			prev_row = cyclone.loc[cyclone['step'] == row['step'] - 1].iloc[0] # localiza o passo de tempo anterior

			d = haversine_np(row['longitude'], row['latitude'], prev_row['longitude'], prev_row['latitude']) 
			distance += d
		
		# atribuindo o resultado às linhas correspondentes
		if distance < 10:
			semi_estacionario.loc[condition] = True

	return semi_estacionario


def haversine_np(lon1, lat1, lon2, lat2):
	"""
	Calculate the great circle distance between two points
	on the earth (specified in decimal degrees)
	
	All args must be of equal length.    
	
	"""
	lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
	
	dlon = lon2 - lon1
	dlat = lat2 - lat1
	
	a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
	
	c = 2 * np.arcsin(np.sqrt(a)) # em radianos
	d = np.rad2deg(c) # em graus
	return d


def grafico_interquartil(data_folder, save_folder):
	df = pd.read_excel(os.path.join(data_folder, "Tracking_2021.xlsx"), header = 0, sheet_name = '2021')

	# Filtrando o que é ciclogenese, o que não é  e o que foi classificado errado.
	df['legenda'] = ""
	df.loc[df['valido'] == 0, 'legenda'] = 'Inválidos'
	df.loc[df['valido'] == 1, 'legenda'] = 'Válidos (Ciclogênese)'
	df.loc[df['valido'].isin([2, 3, 4]), 'legenda'] = 'Válidos (Ciclone maduro ou em oclusão)'

	# removendo campos em branco
	df = df.loc[df['legenda'] != '']

	# filtrando de acordo com a duracao (removendo inferior a 48 horas)
	condition = df['duration'] >= 48
	df = df.loc[condition].copy() # faz uma copia 
	print("Após o filtro de tempo restaram:", np.unique(df['id']).shape[0], "ciclones")

	# filtrando de acordo com o criterio de sistemas semi estacionarios.
	df = df.loc[~sitemas_semiestacionarios(df)].copy()
	
	# Quantos ciclones sobraram depois desses filtros?
	unicos = np.unique(df['id'])
	print(f"Após os dois filtros restaram {unicos.shape[0]} / 402 Ciclones.\nDos quais:\n")
	print(df.groupby(['legenda', 'id']).apply(lambda x: 1).groupby(level = 'legenda').sum())

	# gerando as figuras
	interquartil_vorticidade(df, save_folder)
	interquartil_vorticidade_derivada(df, save_folder)

 