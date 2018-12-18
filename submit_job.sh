#!/bin/bash

### Las líneas #SBATCH configuran los recursos de la tarea
### (aunque parezcan estar comentadas)

### Nombre de la tarea
#SBATCH --job-name=p2c

### Cola de trabajos a la cual enviar.
#SBATCH --partition=batch

### Procesos a largar.
### Es OpenMP, o sea que un proceso en un nodo y muchos hilos.
#SBATCH --ntasks=1

### Hilos por proceso
### Poner el mismo valor acá que en OMP_NUM_THREADS/MKL_NUM_THREADS
#SBATCH --cpus-per-task=56

### Tiempo de ejecucion. Formato dias-horas:minutos.
#SBATCH --time 1-00:00

### Script que se ejecuta al arrancar el trabajo

### Cargar el entorno del usuario incluyendo la funcionalidad de modules
### No tocar
. /etc/profile

### Configurar OpenMP/MKL/etc con la cantidad de cores detectada.
export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=56

### Cargar los módulos para la tarea
# FALTA: Agregar los módulos necesarios
module load clemente

### Largar el programa
# FALTA: Cambiar el nombre del programa
srun python fig3-1.py

