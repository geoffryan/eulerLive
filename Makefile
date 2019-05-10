
BIN = eulerLive

CC = gcc



SDL = /usr/local/Cellar/sdl2/2.0.8

CFLAG = -O3 $(shell sdl2-config --cflags)
SRC = main.c vis.c euler.c
HDR = vis.h euler.h
INC = 
LIB = $(shell sdl2-config --libs) -lm

default: $(BIN)

$(BIN): $(SRC) $(HDR)
	$(CC) $(CFLAG) -o $@ $(SRC) $(INC) $(LIB)

clean:
	rm -f $(BIN)
