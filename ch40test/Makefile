lib = -I./MLDK -L./MLDK -lML64i3  -lpthread -lrt

All: ch40test

pot.tm.c:
	MLDK/mprep src/pot.tm -o src/pot.tm.c
	
pot.tm.o: pot.tm.c
	gcc src/pot.tm.c ${lib} -c

potlib.o:
	gfortran src/potlib.f -c

pot.o:
	gcc src/pot.c ${lib} -c
	
utility.o:
	gfortran src/utility.f -c

ch40test: pot.tm.o potlib.o pot.o utility.o
	gfortran utility.o pot.tm.o potlib.o pot.o ${lib} -o ch40test
	rm *.o
	rm src/*.tm.c
	
clean:
	rm ch40test
	rm *~
