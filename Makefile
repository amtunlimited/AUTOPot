mldk = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
cflags = -I${mldk} -L${mldk} -lML64i3  -lpthread -lrt


All: H3

pot.tm.c:
	${mldk}/mprep src/pot.tm -o src/pot.tm.c
	
pot.tm.o: pot.tm.c
	gcc src/pot.tm.c ${cflags} -c

potlib.o:
	gfortran src/potlib.f -c

pot.o:
	gcc src/pot.c ${cflags} -c
	
utility.o:
	gfortran src/utility.f -c

H3: pot.tm.o potlib.o pot.o utility.o
	gfortran utility.o pot.tm.o potlib.o pot.o ${cflags} -o H3
	rm *.o
	rm src/*.tm.c
	
clean:
	rm H3
	rm *~
