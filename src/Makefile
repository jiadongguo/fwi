CC = mpicc
CFLAGS = -g -Wall
BIN = ../bin
LIB =  -lm -fopenmp -lpthread -lopenblas

LIBS = eal.c fdutils.c cstd.c waveutils.c
INC = -I../include
OBJ = $(SRC:.c=.o)

all: clean fwi smooth ricker
fwi_obj:fwi_obj.c $(LIBS)
	gcc $(CFLAGS) -o $(BIN)/fwi_obj $(INC) $^ $(LIB)
ricker:
	$(CC) $(CFLAGS) -o $(BIN)/ricker $(INC) cstd.c ricker.c $(LIB)
smooth:
	$(CC) $(CFLAGS) -o $(BIN)/smooth $(INC) cstd.c smooth.c $(LIB)
fdmodel:cstd.c eal.c waveutils.c fdutils.c fdmodel.c
	mpicc $(CFLAGS) -o $(BIN)/fdmodel $(INC) cstd.c eal.c waveutils.c fdutils.c fdmodel.c $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) $(BIN)/*
