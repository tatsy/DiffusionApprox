CC		= clang++
SRC		= sources/main.cpp
PROGRAM = diffusion
CFLAGS  = -std=c++11 -O2

all:
	$(CC) $(CFLAGS) -o $(PROGRAM) $(SRC)
