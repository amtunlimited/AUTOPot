mldk = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions

cflags = -I${mldk} -L${mldk} -lML64i3 -lpthread -lrt

ffile = dummy

name = $(notdir $(basename ${ffile}))

All: ${name}

tmp/pot.tm.c:
	cd tmp; \
	${mldk}/mprep pot.tm -o pot.tm.c
	
tmp/pot.tm.o: tmp/pot.tm.c
	cd tmp; \
	gcc pot.tm.c ${cflags} -c -o pot.tm.o

tmp/pot.o:
	gcc src/pot.c ${cflags} -c -o tmp/pot.o
	
tmp/utility.o:
	gfortran src/utility.f -c -o tmp/utility.o
	
tmp/${name}.o:
	gfortran ${ffile} -c -o tmp/${name}.o

${name}: tmp/pot.tm.o tmp/${name}.o tmp/pot.o tmp/utility.o 
	cd tmp; \
	gfortran pot.tm.o ${name}.o utility.o pot.o ${cflags} -o ../${name}
	
clean:
	rm tmp/*
	rm *~ -f
