mldk = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions

cflags = -I${mldk} -L${mldk} -lML64i3 -lpthread -lrt

ffile = dummy

name = $(notdir $(basename ${ffile}))

All: ${name}

pot.tm.c:
	cd tmp; \
	${mldk}/mprep pot.tm -o pot.tm.c
	
pot.tm.o: pot.tm.c
	cd tmp; \
	gcc pot.tm.c ${cflags} -c -o pot.tm.o

${name}.o:
	gfortran ${ffile} -c -o tmp/${name}.o

pot.o:
	gcc src/pot.c ${cflags} -c -o tmp/pot.o
	
utility.o:
	gfortran src/utility.f -c -o tmp/utility.o

${name}: pot.tm.o ${name}.o pot.o utility.o
	cd tmp; \
	gfortran utility.o pot.tm.o ${name}.o pot.o ${cflags} -o ../${name}
	
clean:
	rm tmp/*
	rm *~ -f
