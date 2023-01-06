'''
Arquivo principal para o trabalho de fisica do clima.
'''

# IMPORT BUILT-IN MODULES
import os
import time

# IMPORT CUSTOM MODULES
from Scripts.tracking_postprocessing import post_processing
from Scripts.grafico_derivada import grafico_interquartil
from Scripts import figura_tracking

if __name__ == "__main__":
	folder = os.path.normpath(os.getcwd())
	data_folder = os.path.join(folder, 'Dados')
	figuras_folder = os.path.join(folder, 'Figuras')
	
	# Inicia o timer
	start = time.time()

	# OPERACOES (DEPENDE DO SCRIPT)

	# pos processamento da saida do tracking
	# anos = [2014, 2015, 2017, 2018, 2019, 2021]
	# post_processing(data_folder, anos)

	# Grafico de intervalos interquartis da vorticidade relativa em cada passo de tempo
	xls_folder = os.path.join(data_folder, "Tracking", "Validacao")
	grafico_interquartil(xls_folder, figuras_folder)

	# figura das trajetorias de cada ciclone no ano
	# figura_tracking.leitura_arquivos(folder, 2021)

	# Termina o timer
	end = time.time()

	# Avisa o usuário sobre o término
	print(f"Finalizado!\nTempo de execução: {end - start:.2f}s")

