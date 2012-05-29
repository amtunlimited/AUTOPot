lib = -I./MLDK -L./MLDK -lML32i3 -lpthread -lrt -m32 -mwindows

All: NAME

pot.tm.c:
	MLDK/mprep src/pot.tm -o src/pot.tm.c
	
pot.tm.o: pot.tm.c
	gcc src/pot.tm.c ${lib} -c

potlib.o:
	gfortran src/potlib.f ${lib} -c

pot.o:
	gcc src/pot.c ${lib} -c
	
utility.o:
	gfortran src/utility.f ${lib} -c

NAME: pot.tm.o potlib.o pot.o utility.o
	gfortran utility.o pot.tm.o potlib.o pot.o ${lib} -o NAME
	rm *.o
	rm src/*.tm.c
	
clean:
	rm NAME
	rm *~
