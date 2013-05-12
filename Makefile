

#For debugging:
#CFLAGS=-Wall -Werror -g

#For actually running:
FLAGS=-Wall -Werror -O2


all: 2dfold 3dfold 

2dfold: proteins.o list.o

3dfold: proteins.o list.o

clean:
	rm -f *.o *~ test 2dfold 3dfold output
