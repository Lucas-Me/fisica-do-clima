'''
Arquivo principal para o trabalho de fisica do clima.
'''

# IMPORT BUILT-IN MODULES
import os
import time

# IMPORT MODULES
import numpy as np

# IMPORT CUSTOM MODULES
from Scripts.tracking_postprocessing import post_processing
from Scripts.grafico_derivada import grafico_interquartil
from Scripts import figura_tracking
from Scripts import densidade_ciclogenetica
from Scripts.graficos_espaciais import sazonalidade_ciclogenese, visualizar_mapas
from Scripts.serie_temporal import ciclogenese_temporal
from Scripts.mapas_era5 import tendencia_fluxos

if __name__ == "__main__":
	# diretorios
	folder = os.path.normpath(os.getcwd())
	data_folder = os.path.join(folder, 'Dados')
	figuras_folder = os.path.join(folder, 'Figuras')
	dados_filtrados_dir = os.path.join(data_folder, 'Tracking', 'Filtrado')
	processados_dir = os.path.join(data_folder, 'Tracking', 'Processado')
	#
	era5_dir = os.path.join(data_folder, 'ERA5')
	
	# configuracoes
	anos = np.linspace(1991, 2021, 31).astype(int)

	# Inicia o timer
	start = time.time()

	# POS PROCESSAMENTO
	# ///////////////////////////////////////////////////////////////////////

	# pos processamento da saida do tracking
	# post_processing(data_folder, anos)

	# Grafico de intervalos interquartis da vorticidade relativa em cada passo de tempo
	# xls_folder = os.path.join(data_folder, "Tracking", "Validacao")
	# grafico_interquartil(xls_folder, figuras_folder)

	# figura das trajetorias de cada ciclone no ano (MT FIGURA, ESTEJA AVISADO)
	# figura_tracking.leitura_arquivos(folder, 2021)

	# DENSIDADE CICLOGENETICA
	# /////////////////////////////////////////////////////////////////////

	# Confecção do NetCDF com a densidade ciclogenética
	# densidade_ciclogenetica.gerar_netcdf(dados_filtrados_dir, anos, resolucao = 1)

	# Visualizar plot basico do netcdf gerado
	# visualizar_mapas(processados_dir, figuras_folder)

	# # serie temporal da quantidade de ciclogenese por regiao ciclogenetica
	# densidade_ciclogenetica.regioes_ciclogeneticas(dados_filtrados_dir, anos)
	# densidade_ciclogenetica.areas_ciclogeneticas(dados_filtrados_dir, anos, sazonal = True)
	# densidade_ciclogenetica.areas_ciclogeneticas(dados_filtrados_dir, anos)

	# # Sazonalidade espacial ciclogenese
	# sazonalidade_ciclogenese(processados_dir, figuras_folder)

	# SERIES TEMPORAIS
	# ///////////////////////////////////////////////////////////////////////
	# ciclogenese_temporal(processados_dir, figuras_folder)

	# FLUXOS
	# //////////////////////////////////////////////////////////////////////
	tendencia_fluxos(era5_dir)
	
	# Termina o timer
	end = time.time()

	# Avisa o usuário sobre o término
	print(f"Finalizado!\nTempo de execução: {end - start:.2f}s")

