

CFLAGS =  -std=c++11 -g
#CFLAGS =  -g
CC = g++
AR = ar
ARFLAGS = rv
MYLIB = libmine
LIB =
LINK =   -lm
.SUFFIXES: .c .cpp



.cpp.o:
	echo 'Compiling $<'
	$(CC) $(CFLAGS) -c $<


.c.o:
	$(CC) $(CFLAGS) -c $<


all:	systemSolver 


systemSolver:	systemSolver.o poisson.h poisson.o solution.o solution.h preconditioner.o preconditioner.h 
	echo $@
	$(CC) -o $@ $@.o  poisson.o solution.o preconditioner.o $(LINK) 


clean:	
	rm -f *.o systemSolver 



