lib = -I./MLDK -L./MLDK -lml64i3 -m64

All: NAME

pot.tm.c:
	MLDK/mprep src/pot.tm -o src/pottm.c
	
pottm.o: pot.tm.c
	gcc src/pottm.c ${lib} -c

potlib.o:
	gfortran src/potlib.f ${lib} -c

pot.o:
	gcc src/pot.c ${lib} -c
	
utility.o:
	gfortran src/utility.f ${lib} -c

NAME: pottm.o potlib.o pot.o utility.o
	gfortran utility.o potlib.o pot.o pottm.o ${lib} -o NAME
	rm *.o
	rm src/*.tm.c
	
clean:
	rm NAME
	rm *~
