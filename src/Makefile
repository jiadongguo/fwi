CC = mpicc
CFLAGS = -g -Wall -I../include
BIN = ../bin
LIB =  -lm -fopenmp -lpthread -lopenblas
LIBS = cstd.c eal.c lap.c waveutils.c
all: clean fdmodeling2 mkmodel fwi smooth ricker

%:%.c $(LIBS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(LIB)

# fdmodeling2:$(OBJ) main.c
# 	$(CC) $(CFLAGS) -o $(BIN)/fdmodeling2 $(OBJ) $(LIB)
clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -rf $(BIN)
	mkdir $(BIN)