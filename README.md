# metodos_numericos
Programas de la materia de Mario, métodos numéricos en astrofisica

para compilar el .f90 tenemos que poner 
f2py -c -m aceleracion aceleracion.f90

f2py --fcompiler=gfortran --opt='-fopenmp -O3' -c -m integradorespdc integradores_modificado.f90


Para mulatona:
f2py --fcompiler=gfortran --opt='-fopenmp -O3 -fPIC -march=broadwell' -c -m ipdc integ_cap3.f90