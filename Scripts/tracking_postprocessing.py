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
	timedelta = pd.Timedelta(1, 'hour')
	dt = pd.Timedelta(6, 'hour') # duracao de um passo de tempo

	# loop para cada id / ciclone
	for n in ids:
		# slice do dataframe para o tracking do ciclone em questao
		condition = tracking['id'] == n
		cyclone_dt = tracking.loc[condition, 'datetime'].copy()

		# extrai a data mais recente e mais atiga
		min_dt = cyclone_dt.min()
		max_dt = cyclone_dt.max()

		# estima a duracao em hora
		duracao = (max_dt - min_dt) / timedelta

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
	'''

	# id's unicos
	tracking = tracking.copy()
	ids = np.unique(tracking['id'])
	
	# timedelta de cada passo de tempo
	dt_h = 6 

	# derivadas
	derivada_1 = pd.Series('nan', index = tracking.index)
	derivada_2 = pd.Series('nan', index = tracking.index)

	# loop para cada id / ciclone
	for n in ids:
		# slice do dataframe para o tracking do ciclone em questao
		condition = tracking['id'] == n
		cyclone = tracking.loc[condition].copy()

		# organizando em ordem crescente de acordo com o tempo
		cyclone = cyclone.sort_values(by = 'datetime', ascending = True)
		vort = cyclone['rel_vort']

		# derivada temporal
		dvort_dt = [(vort.iloc[i + 1] - vort.iloc[i]) / dt_h for i in range(vort.shape[0] - 1)] 
		dvort_dt2 = [(dvort_dt[i + 1] - dvort_dt[i]) / dt_h for i in range(len(dvort_dt) - 1)] 

		derivada_1.loc[cyclone.index] = dvort_dt + ['nan']
		derivada_2.loc[cyclone.index] = dvort_dt2 + ['nan'] * 2

	tracking['derivada_1'] = derivada_1
	tracking['derivada_2'] = derivada_2
	
	return tracking


def post_processing(folder : str, years : list[int]) -> None:
	# Diretorio dos dados de tracking
	data_folder = os.path.join(folder, "Tracking")

	# criando diretorio para guardar os dados processados, se ja nao existe.
	save_folder = os.path.join(data_folder, "Processado")
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

		# estimando a duracao de cada ciclone
		# ////////////////////////////////////////////////////////////////////
		df = cyclone_duration(df)

		# Quantos obedecem ao criterio?
		# n = df.shape[0] - df['criterio_dt'].sum() 
		# print(f"Em {year}, houveram {n} ({n / df.shape[0] * 100:.2f}%) ciclones com duração inferior a 24 horas")

		# estimando o criterio de vorticidade relativa
		# ////////////////////////////////////////////////////////////////////
		df = relative_vorticity(df)

		# Salvando
		save_file = os.path.join(save_folder, f"{year}.csv")
		df.to_csv(save_file, index= False, header = True)