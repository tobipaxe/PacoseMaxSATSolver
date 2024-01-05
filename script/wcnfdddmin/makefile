CC = g++
LIBS = -lz
CFLAGS = -O3 -DNDEBUG -std=c++17
DEBUGCFLAGS = -g -Wall -Wextra -DLOGGING
NAME = wcnfddmin

all:
	$(CC) $(CFLAGS) $(NAME).cpp -o $(NAME) $(LIBS)

debug:
	$(CC) $(DEBUGCFLAGS) $(NAME).cpp -o $(NAME) $(LIBS)

clean:
	rm $(NAME)
