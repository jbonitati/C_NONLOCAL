#
CC = gcc

#edit LAPACK_PATH if necessary
LAPACK_PATH = /usr/lib64/atlas

a,out: main.c
$(CC)  main.c  -L$(LAPACK_PATH) -llapack -lblas  -lgfortran -lm


