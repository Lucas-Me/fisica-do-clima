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

def cyclone_duration(tracking : pd.DataFrame):
	'''
	Estima a duração de cada ciclone rastreado.
	Para isso, recebe um DataFrame com os dados em tabela, realiza um loop no id
	de cada ciclone e retorna uma nova coluna com a duracao em horas.	
	'''

	# id's unicos
	tracking = tracking.copy()
	ids = np.unique(tracking['id'])
	
	# cria uma serie com valores temporarios
	duration = pd.Series(0, index = tracking.index)
	passo = pd.Series(0, index = tracking.index)
	
	# timedelta
	dt = pd.Timedelta(6, 'hour') # duracao de um passo de tempo

	# loop para cada id / ciclone
	for n in ids:
		# slice do dataframe para o tracking do ciclone em questao
		condition = tracking['id'] == n
		cyclone_dt = tracking.loc[condition, 'datetime'].copy()

		# extrai a data mais recente e mais atiga
		min_dt = cyclone_dt.min()
		max_dt = cyclone_dt.max()

		# estima a duracao em horas
		duracao = (max_dt - min_dt) / pd.Timedelta(1, 'hour')

		# atualiza a série
		duration.loc[condition] = duracao
		passo.loc[condition] = (cyclone_dt - min_dt) // dt

	tracking['duration'] = duration
	tracking['step'] = passo

	return tracking


def relative_vorticity(tracking: pd.DataFrame):
	'''
	Estima o criterio de vorticidade relativa para cada ciclone
	Para isso, recebe um DataFrame com os dados em tabela, realiza um loop no id
	de cada ciclone e retorna uma nova coluna com o resultado (True ou False).

	O criterio define que o ciclone sera filtrado se após 72h do 1º passo de tempo
	em que foi detectado pelo tracking, o seu valor de vorticidade relativa (-)
	for inferior a 4.

	Estima-se que isso elimine 90% das classificacoes erradas, 50% dos ciclones
	ja em estado maduro / oclusao e 25% das ciclogeneses.

	Retorna True se for inferior a 4, False caso contrario.
	'''

	# id's unicos
	tracking = tracking.copy()
	ids = np.unique(tracking['id'])
	
	# derivadas
	criterio = pd.Series(False, index = tracking.index)

	# passo de tempo
	dt_h = 6 # resolucao temporal em horas
	dt_criterio = 60
	step = dt_criterio // dt_h

	# loop para cada id / ciclone
	for n in ids:
		# slice do dataframe para o tracking do ciclone em questao
		condition = tracking['id'] == n
		cyclone = tracking.loc[condition]

		if cyclone['duration'].iloc[0] < dt_criterio: # se a duracao do ciclone for inferior a 72h, pula.
			continue

		# organizando em ordem crescente de acordo com o tempo
		vort = cyclone.loc[cyclone['step'] == step, 'rel_vort'].iloc[0]
		# ini_vort = cyclone.loc[cyclone['step'] == 0, 'rel_vort'].iloc[0]

		if vort < 4:
		# if ini_vort < 2.5:
			criterio.loc[condition] = True
	
	return criterio


def post_processing(folder : str, years : list[int]) -> None:
	# Diretorio dos dados de tracking
	data_folder = os.path.join(folder, "Tracking")

	# criando diretorio para guardar os dados processados, se ja nao existe.
	save_folder = os.path.join(data_folder, "Filtrado")
	os.makedirs(save_folder, exist_ok = True)

	# Colunas do arquivo
	cols = ['id', 'datetime', 'longitude', 'latitude', 'rel_vort', 'vazia']
	
	# tipos de dados de cada coluna
	types = [int, str, float, float, float, str]
	dtypes = dict(zip(cols, types))

	# argumentos da funcao read_csv
	kwargs = dict(
		names = cols,
		dtype = dtypes,
		decimal = '.',
		sep = ',',
		index_col = False,
		header = None
	)

	# Lendo os arquivos
	for year in years:
		file_path = os.path.join(data_folder, f'{year}.csv')
		df = pd.read_csv(file_path, **kwargs)

		# dropando a coluna vazia
		df = df.drop(columns  = 'vazia')

		# convertendo datetime
		df['datetime'] = df['datetime'].str.slice(2, -1)
		df['datetime'] = pd.to_datetime(df['datetime'])

		N = np.unique(df['id']).shape[0] # quantidade original de dados

		# estimando a duracao de cada ciclone e criterio de 48h
		# ////////////////////////////////////////////////////////////////////
		df = cyclone_duration(df)
		df = df.loc[df['duration'] >= 48]

		# estimando o criterio de vorticidade relativa
		# ////////////////////////////////////////////////////////////////////
		df = df.loc[~relative_vorticity(df)] # True se rel_vort for inferior a 4 em 60h


		print(f"Após os filtros, restaram {np.unique(df['id']).shape[0]} / {N} ciclones em {year}")
		# Salvando
		save_file = os.path.join(save_folder, f"{year}.csv")
		df.to_csv(save_file, index= False, header = True)