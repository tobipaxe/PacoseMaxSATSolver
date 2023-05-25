CC = g++
CFLAGS = -std=c++11 -O2 -Wall -Wextra -Wshadow -g -Wfatal-errors -I.
OBJFILES = utils/System.o
LIBFILES =


.PHONY: lib

lib: core/Solver.o $(OBJFILES)
	mkdir -p ../lib
	ar rvs ../lib/libglucose3.a core/Solver.o $(OBJFILES) $(SATOBJFILES)

%.o: %.cc
	$(CC) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -f ../lib/libglucose3.a
	rm -f *.o *.d
	rm -f core/*.o core/*.d
	rm -f utils/*.o utils/*.d
