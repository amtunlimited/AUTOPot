mldk = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions

cflags = -I${mldk} -L${mldk} -lML64i3 -lpthread -lrt

file = dummy

All: ${file}

pot.tm.c:
	${mldk}/mprep src/pot.tm -o src/pot.tm.c
	
pot.tm.o: pot.tm.c
	gcc src/pot.tm.c ${cflags} -c

${file}.o:
	gfortran src/${file}.f -c

pot.o:
	gcc src/pot.c ${cflags} -c
	
utility.o:
	gfortran src/utility.f -c

${file}: pot.tm.o ${file}.o pot.o utility.o
	gfortran utility.o pot.tm.o ${file}.o pot.o ${cflags} -o ${file}
	rm *.o
	rm src/*.tm.c
	
clean:
	rm H3
	rm *~
