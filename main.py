'''
Arquivo principal para o trabalho de fisica do clima.
'''

# IMPORT BUILT-IN MODULES
import os
import time

# IMPORT CUSTOM MODULES
from Scripts import 


if __name__ == "__main__":
	folder = os.path.normpath(os.getcwd())
	data_folder = os.path.join(folder, 'Dados')
	
	# Inicia o timer
	start = time.time()

	# Operacoes (Depende do script)


	# Termina o timer
	end = time.time()


	# Avisa o usuário sobre o término
	print(f"Finalizado!\nTempo de execução: {end - start:.2f}")

