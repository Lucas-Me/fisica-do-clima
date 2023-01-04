'''
Arquivo principal para o trabalho de fisica do clima.
'''

# IMPORT BUILT-IN MODULES
import os
import time

# IMPORT CUSTOM MODULES
from Scripts.tracking_postprocessing import post_processing
from Scripts.grafico_derivada import grafico_interquartil

if __name__ == "__main__":
	folder = os.path.normpath(os.getcwd())
	data_folder = os.path.join(folder, 'Dados')
	
	# Inicia o timer
	start = time.time()

	# OPERACOES (DEPENDE DO SCRIPT)

	# pos processamento da saida do tracking
	# anos = [2014, 2015, 2017, 2018, 2019, 2021]
	# post_processing(data_folder, anos)

	# Grafico de intervalos interquartis da derivada temporal em cada passo de tempo
	grafico_interquartil(folder, 2014)

	# Termina o timer
	end = time.time()

	# Avisa o usuário sobre o término
	print(f"Finalizado!\nTempo de execução: {end - start:.2f}s")

