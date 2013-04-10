CFLAGS=-Wall -Werror


all: 2dfold 3dfold test

2dfold: proteins.o list.o

3dfold: proteins.o list.o

test: proteins.o list.o

clean:
	rm -f *.o test 2dfold 3dfold output
