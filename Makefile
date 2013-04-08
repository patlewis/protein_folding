CFLAGS=-Wall -Werror

2dfold: proteins.o 2d_folding.c
	cc -O2 proteins.c 2d_folding.c

3dfold: proteins.o 3d_folding.c
	cc -O2 proteins.c 3d_folding.c

proteins.o:
	cc -c proteins.c

all:
	2d
	3d
