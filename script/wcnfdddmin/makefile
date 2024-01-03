CC = g++
LIBS = -lz
CFLAGS = -O3 -DNDEBUG
DEBUGCFLAGS = -g -Wall -Wextra -DLOGGING
NAME = wcnfddmin

all:
	$(CC) $(CFLAGS) $(NAME).cpp -o $(NAME) $(LIBS)

debug:
	$(CC) $(DEBUGCFLAGS) $(NAME).cpp -o $(NAME) $(LIBS)

clean:
	rm $(NAME)
