cflags = -I/usr/include -L/home/arn44621/.Wolfram/9.0/SystemFiles/Links/MathLink/DeveloperKit/Linux/CompilerAdditions -I/home/arn44621/.Wolfram/9.0/SystemFiles/Links/MathLink/DeveloperKit/Linux/CompilerAdditions -lML64i3 -lpthread -lrt

ffile = fh25sec.f

name = $(notdir $(basename ${ffile}))

All: ${name}

bin/pot.tm.o:
	gcc src/pot.tm.c ${cflags} -c -o bin/pot.tm.o

bin/pot.o:
	gcc src/pot.c ${cflags} -c -o bin/pot.o
	
bin/utility.o:
	gfortran src/utility.f -c -o bin/utility.o
	
bin/${name}.o:
	gfortran src/${ffile} -c -o bin/${name}.o

${name}: bin/pot.tm.o bin/${name}.o bin/pot.o bin/utility.o 
	cd bin; \
	gfortran pot.tm.o ${name}.o utility.o pot.o ${cflags} -o ../${name}
	
clean:
	rm src/*
	rm *~ -f
