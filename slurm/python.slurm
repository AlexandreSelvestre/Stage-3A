#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --partition=cpu_med
#SBATCH --job-name=extrac_multi
#SBATCH --output=%x.o%j
#SBATCH --export=PATH,LD_LIBRARY_PATH,FONTCONFIG_PATH,HOME,USER

echo "[$(date)] Début du script"

module purge  
echo "[$(date)] Modules purgés"

module load anaconda3/2024.06/gcc-13.2.0
echo "[$(date)] Module anaconda3 chargé"

source activate radio
echo "[$(date)] Environnement conda 'radio' activé"

cd /gpfs/users/selvestra/python_codes
echo "[$(date)] Répertoire changé"

echo "[$(date)] Début de l'exécution du script Python"
python utils_advanced_area.py
echo "[$(date)] Fin de l'exécution du script Python"
