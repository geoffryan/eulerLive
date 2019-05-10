
BIN = eulerLive

CC = gcc

SDL = /usr/local/Cellar/sdl2/2.0.8

CFLAG = -O3
SRC = main.c vis.c euler.c
HDR = vis.h euler.h
INC = -I$(SDL)/include
LIB = -L$(SDL)/lib -lSDL2 -lm

default: $(BIN)

$(BIN): $(SRC) $(HDR)
	$(CC) $(CFLAG) -o $@ $(SRC) $(INC) $(LIB)

clean:
	rm -f $(BIN)
