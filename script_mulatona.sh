#!/bin/bash

### Nombre de la tarea
#SBATCH --job-name=p2c
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time 1-00:00
### Cargar el entorno del usuario incluyendo la funcionalidad de modules

### No tocar
. /etc/profile

### Configurar OpenMP/MKL/etc con la cantidad de cores detectada.
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export OMP_STACKSIZE=5000M

### Cargar los módulos para la tarea
module load  gcc/7
srun python fig3-1.py
