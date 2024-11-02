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
ricker:ricker.c ${LIBS}
	$(CC) $(CFLAGS) -o $(BIN)/ricker $(INC) $^ $(LIB)
smooth:
	$(CC) $(CFLAGS) -o $(BIN)/smooth $(INC) cstd.c smooth.c $(LIB)
fdmodel:fdmodel.c ${LIBS}
	mpicc $(CFLAGS) -o $(BIN)/fdmodel $(INC) $^ $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) $(BIN)/*
