GITV = $(shell git describe --always --dirty)
CC = g++
CFLAGS = -std=c++11 -O2 -Wall -Wextra -Wshadow -g -Wfatal-errors -fPIC -Isatsolver/solvers/glucose-3.0/ -DGIT_IDENTIFIER=\"$(GITV)\"
OBJFILES = preprocessor.o inputreader.o outputreader.o preprocessedinstance.o trace.o utility.o probleminstance.o timer.o clause.o log.o AMSLEX.o touchedlist.o preprocessorinterface.o cardinalityconstraint.o satlikeinterface.o
LIBFILES = satsolver/solvers/lib/libglucose3.a
SATOBJFILES = satsolver/solvers/glucose3/core/Solver.o satsolver/solvers/glucose3/utils/System.o

all: preprocessor

preprocessor: satsolvers main.o $(OBJFILES)
	$(CC) $(CFLAGS) main.o $(OBJFILES) $(LIBFILES) -o maxpre

.PHONY: lib
lib: satsolvers preprocessorinterface.o $(OBJFILES)
	mkdir -p lib
	ar rvs lib/libmaxpre.a preprocessorinterface.o $(OBJFILES) $(SATOBJFILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -MMD -c $< -o $@

-include *.d

satsolvers: glucose3

glucose3:
		$(MAKE) -C satsolver/solvers/glucose3 lib

clean:
	rm -f lib/*.a
	rm -f *.o *.d
	$(MAKE) -C satsolver/solvers/glucose3 clean
