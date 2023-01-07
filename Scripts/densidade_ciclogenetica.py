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
	lon_max = 0
	lat_min = -70
	lat_max = 20

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

	# Divisao pela area corrspondente a cada pixel. 
	# Calculo de um elemento de area dA no globo da Terra
	# Terra foi considerada uma esfera perfeita (nao é :p)
	r = 6371 # raio da terra em kilometros
	dlon = dx * np.pi / 180 # intervalo de longitude
	lats = Y - dx / 2 # latitudes
	dA = np.abs(r ** 2 * dlon * ( np.sin( np.deg2rad( lats + dx ) ) - np.sin( np.deg2rad( lats ) ) ) ) # em km²

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

		#  indica a  data inicial do ciclone em cada passo de tempo
		df['mes_inicio'] = pd.Series(0, df.index)
		ids_unicos = np.unique(df['id'])
		for id_ in ids_unicos:
			condition = df['id'] == id_
			cyclone = df.loc[condition]
			mes = cyclone['datetime'].dt.month.iloc[0]
			df.loc[condition, 'mes_inicio'] = mes

		# loop para cada mes
		for month in meses:
			# slice para o mes e ano
			subset = df.loc[df['mes_inicio'] == month]

			# resgatando a media dos primeiros passos de tempo
			group = subset.loc[subset['step'] <= 4].groupby('id')
			grouped = group.agg(
				lon = ('longitude', 'mean'),
			)
			lat = subset.loc[subset['step'] == 0, ['latitude', 'id']]

			# calculando a densidade
			idx = (year - ano_inicio) * 12 + (month - 1)
			results[:, :, idx] =  matriz_densidade(
				binx,
				biny,
				grouped['lon'].to_numpy(),
				lat['latitude'].to_numpy()
				) / dA # Dividindo pelos infinitesimos de area

	# Gerando o arquivo netcdf
	da = xr.DataArray(results, coords = coords, dims = dims)
	ds = xr.Dataset({"ciclogenese" : da})
	filename = 'frequencia_ciclogenetica.nc'
	ds.to_netcdf(os.path.join(save_folder, filename))


def areas_ciclogeneticas(data_folder, years, sazonal = False):
	'''
	Gera um arquivo csv com a série temporal da quantidade de ciclogeneses
	por areas ciclogenetica, conforme defeinida através da analise da figura
	de densidade de ciclogeneses.

	Input são os ciclones identificdos pelo tracking após o pos-processamento 
	(aplicacao do filtro)

	Years deve ser uma lista em ordem crescente.
	'''

	# cria uma pasta "Processado" no diretorio pai, tudo bem se ja houver
	par_dir = os.path.abspath(os.path.join(data_folder, os.pardir))
	save_folder = os.path.join(par_dir, 'Processado')
	os.makedirs(save_folder, exist_ok=True)

	# array de datas
	times = pd.date_range(start = f"{years[0]}-1-1", end = f"{years[-1]}-12-1", freq = 'MS') 

	#A1 (Regiao Ciclogenetica no litoral do Urugai/ regiao sul)
	A1_lon = slice(-60.0, -50.0)
	A1_lat = slice(-27.5, -37.5)
	A1_values = np.zeros(times.shape)

	#A2 (Regiao ciclogenetica no litoral da Argentina)
	A2_lon = slice(-62.5, -52.5)
	A2_lat = slice(-40.0, -50.0)
	A2_values = np.zeros(times.shape)

	#A3 (Regiao ciclogenetica na patagonia)
	A3_lon = slice(-72.5, -62.5)
	A3_lat = slice(-45.0, -55.0)
	A3_values = np.zeros(times.shape)

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

		# estimando a qual regiao ciclogeneticapertence cada ponto
		df['regiao'] = "Nenhuma"

		df.loc[((df['longitude'] >= A1_lon.start) & (df['longitude'] <= A1_lon.stop))
							& ((df['latitude'] >= A1_lat.stop) & (df['latitude'] <= A1_lat.start)), 'regiao'] = "A1"

		df.loc[((df['longitude'] >= A2_lon.start) & (df['longitude'] <= A2_lon.stop))
							& ((df['latitude'] >= A2_lat.stop) & (df['latitude'] <= A2_lat.start)), 'regiao'] = "A2"
			
		df.loc[((df['longitude'] >= A3_lon.start) & (df['longitude'] <= A3_lon.stop))
							& ((df['latitude'] >= A3_lat.stop) & (df['latitude'] <= A3_lat.start)), 'regiao'] = "A3"

		# loop para cada mes
		for month in meses:
			# slice para o mes e ano
			subset = df.loc[df['datetime'].dt.month == month]

			# resgatando o primeiro passo de tempo de cada ciclone
			subset = subset.loc[subset['step'] == 0]

			# estimando a quantidade de pontos (ciclogeneses) por regiao
			qtd = subset['regiao'].groupby(subset['regiao']).size()
			qtd = qtd.reindex(['Nenhuma', 'A1', 'A2', 'A3'], fill_value= 0)

			# posicao no array
			idx = (year - ano_inicio) * 12 + (month - 1) 

			# atribuindo as quantidades
			A1_values[idx] = qtd['A1']
			A2_values[idx] = qtd['A2']
			A3_values[idx] = qtd['A3']

	# Criando um novo DataFrame para armazenar os resultados
	content = {'time' : times, 'A1' : A1_values, 'A2' : A2_values, 'A3' : A3_values}
	new_df = pd.DataFrame.from_dict(content)
	new_df['ano'] = new_df.time.dt.year
	new_df['month'] = new_df.time.dt.month

	# salvando
	filename = 'frequencia_ciclogenese_mensal_AREAS.csv'
	if sazonal:
		new_df['trimestre'] = (new_df['month'] // 3) % 4
		new_df['estacao'] = new_df['trimestre'].apply(lambda x: ['VERAO', 'OUTONO', 'INVERNO', 'PRIMAVERA'][x])
		new_df['trimestre'] = new_df['trimestre'].apply(lambda x: ['DJF', 'MAM', 'JJA', 'SON'][x])
		new_df = new_df.drop(['month', 'time'], axis = 1).groupby(['trimestre', 'ano']).agg(
			A1 = ('A1', 'sum'),
			A2 = ('A2', 'sum'),
			A3 = ('A3', 'sum'),
		)

		new_df = new_df.reset_index()
		filename = 'frequencia_cilcogenese_sazonal_AREAS.csv'

	new_df.to_csv(os.path.join(save_folder, filename), index = False)
	

def regioes_ciclogeneticas(data_folder, years):
	'''
	Gera um arquivo csv com a série temporal da quantidade de ciclogeneses
	por regiao ciclogenetica, conforme definido por Reboita et a., (20xx).

	Input são os ciclones identificdos pelo tracking após o pos-processamento 
	(aplicacao do filtro)

	Years deve ser uma lista em ordem crescente.
	'''

	# cria uma pasta "Processado" no diretorio pai, tudo bem se ja houver
	par_dir = os.path.abspath(os.path.join(data_folder, os.pardir))
	save_folder = os.path.join(par_dir, 'Processado')
	os.makedirs(save_folder, exist_ok=True)

	# array de datas
	times = pd.date_range(start = f"{years[0]}-1-1", end = f"{years[-1]}-12-1", freq = 'MS') 

	#RG1 (Regiao Ciclogenetica no litoral do Brasil sul/sudeste)
	RG1_lon = slice(-49., -36.)
	RG1_lat = slice(-24., -33.)
	RG1_values = np.zeros(times.shape)

	#RG2 (Regiao ciclogenetica no litoral do uruguai)
	RG2_lon = slice(-59., -46.)
	RG2_lat = slice(-34., -43.)
	RG2_values = np.zeros(times.shape)

	#RG3 (Regiao ciclogenetica no litoral da argentina)
	RG3_lon = slice(-68., -55.)
	RG3_lat = slice(-43., -52.)
	RG3_values = np.zeros(times.shape)

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

		# estimando a qual regiao ciclogeneticapertence cada ponto
		df['regiao'] = "Nenhuma"

		df.loc[((df['longitude'] >= RG1_lon.start) & (df['longitude'] <= RG1_lon.stop))
							& ((df['latitude'] >= RG1_lat.stop) & (df['latitude'] <= RG1_lat.start)), 'regiao'] = "RG1"

		df.loc[((df['longitude'] >= RG2_lon.start) & (df['longitude'] <= RG2_lon.stop))
							& ((df['latitude'] >= RG2_lat.stop) & (df['latitude'] <= RG2_lat.start)), 'regiao'] = "RG2"
			
		df.loc[((df['longitude'] >= RG3_lon.start) & (df['longitude'] <= RG3_lon.stop))
							& ((df['latitude'] >= RG3_lat.stop) & (df['latitude'] <= RG3_lat.start)), 'regiao'] = "RG3"

		# loop para cada mes
		for month in meses:
			# slice para o mes e ano
			subset = df.loc[df['datetime'].dt.month == month]

			# resgatando o primeiro passo de tempo de cada ciclone
			subset = subset.loc[subset['step'] == 0]

			# estimando a quantidade de pontos (ciclogeneses) por regiao
			qtd = subset['regiao'].groupby(subset['regiao']).size()
			qtd = qtd.reindex(['Nenhuma', 'RG1', 'RG2', 'RG3'], fill_value= 0)

			# posicao no array
			idx = (year - ano_inicio) * 12 + (month - 1) 

			# atribuindo as quantidades
			RG1_values[idx] = qtd['RG1']
			RG2_values[idx] = qtd['RG2']
			RG3_values[idx] = qtd['RG3']

	# Criando um novo DataFrame para armazenar os resultados
	content = {'time' : times, 'l1' : RG1_values, 'l2' : RG2_values, 'l3' : RG3_values}
	new_df = pd.DataFrame.from_dict(content)
	new_df['ano'] = new_df.time.dt.year
	new_df['month'] = new_df.time.dt.month

	# salvando
	filename = 'frequencia_ciclogenese_mensal.csv'
	new_df.to_csv(os.path.join(save_folder, filename), index = False)
	