"""
Este arquivo trata-se de um dos scripts criados para a confecção do seminário de Física do Clima.

Possui como objetivo criar uma figura com os intervalos interquartis da derivada temporal em cada passo de tempo dos ciclones.

Criado por: Lucas Menezes
"""

# IMPORT MODULES
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def grafico_interquartil(folder, year):
	data_folder = os.path.join(folder, "Dados\\Tracking\\Processado")
	save_folder = os.path.join(folder, "Figuras")

	# Lendo o arquivo processado
	df = pd.read_csv(os.path.join(data_folder, f"{year}.csv"), header = 0)

	# Intervalos interquartis nas primeiras 48h
	group = df.groupby(["step"]) # agrupando por passo de tempo
	
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
	fig, ax = plt.subplots(figsize = (10, 6))

	# Plotando os resultados

	# Mediana  e media
	ax.plot(results.index, results['P50'], color = 'green', linewidth = 2)
	ax.plot(results.index, results['MEDIA'], color = 'royalblue', linewidth = 2)

	# preenche entre o percentil 25 e 75
	ax.fill_between(results.index, results['P25'], results['P75'], color = 'orangered')

	# preenche entre o percentil 5 e 25
	ax.fill_between(results.index, results['P5'], results['P25'], color = 'orange')

	# preenche entre o percentil 75 e 95
	ax.fill_between(results.index, results['P75'], results['P95'], color = 'orange')

	# FFORMATACAO DO EIXO X.
	ax.set_xticks(results.index)
	ax.set_xticklabels(["t"] + [f"t + {passo}" for passo in results.index[1:]])
	ax.set_xlim(0, 7) # Primeiras 48 horas
	ax.set_xlabel("Passo de tempo", fontsize = 12)

	# Formatacao do EIXO Y.
	ax.set_ylabel("Derivada Temporal da Vorticidade Relativa", fontsize = 12)

	# Titulo
	ax.set_title("Ciclo de vida dos ciclones (2014)", fontsize = 18)

	# Demais configuracoes
	ax.grid(True, axis = 'y')

	# Salva a figura e fecha
	nome = os.path.join(save_folder, "Intervalos Interquartis.png")
	print(nome)
	fig.savefig(nome, bbox_inches = 'tight', dpi = 200)
	plt.close()
