CC = mpicc
CFLAGS = -g -Wall
BIN = ../bin
LIB =  -lm -fopenmp -lpthread -lopenblas
INC = -I../include
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: clean fdmodeling2

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

fdmodeling2:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/fdmodeling2 $(OBJ) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) $(BIN)/test
