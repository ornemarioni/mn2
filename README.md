# metodos_numericos
Programas de la materia de Mario, métodos numéricos en astrofisica

para compilar el .f90 tenemos que poner 
f2py -c -m aceleracion aceleracion.f90
f2py --fcompiler=gfortran --opt='-fopenmp -O3' -c -m integradores2 integradores2.f90
