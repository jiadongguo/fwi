CC = mpicc
CFLAGS = -g -Wall
BIN = ../bin
LIB =  -lm -fopenmp -lpthread -lopenblas
INC = -I../include
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: clean fwi smooth

smooth:
	$(CC) $(CFLAGS) -o $(BIN)/smooth $(INC) cstd.c smooth.c $(LIB)
fwi:
	$(CC) $(CFLAGS) -o $(BIN)/fwi $(INC) cstd.c eal.c waveutils.c main.c $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) $(BIN)/*