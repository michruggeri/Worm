LIB  = Random.o ran2.o Link.o Potentials.o System.o worm_canonical.o
CC   = g++
OPT  = -O3 -w
MAIN = Worm

All: ${MAIN}.x 

${MAIN}.x: ${LIB}
	${CC} ${OPT} ${LIB} -o ${MAIN}.x

ran2.o: ran2.hpp ran2.cpp
	${CC} ${OPT} ran2.cpp -c

Random.o: Random.hpp Random.cpp ran2.o
	${CC} ${OPT} Random.cpp -c

Link.o: Link.hpp Link.cpp
	${CC} ${OPT} Link.cpp -c

Potentials.o: Potentials.hpp Potentials.cpp
	${CC} ${OPT} Potentials.cpp -c

System.o: System.hpp System.cpp Potentials.o Link.o Random.o
	${CC} ${OPT} System.cpp -c

worm_canonical.o: worm_canonical.cpp System.o
	${CC} ${OPT} worm_canonical.cpp -c

clean:
	rm -f *.o *.x *.tgz

backup:
	rm -rf ./BCKP; mkdir ./BCKP; cp *.cpp *.hpp Makefile input.example ./BCKP

zip:
	rm worm.tgz; tar -cvzf worm.tgz ran2.cpp ran2.hpp Random.cpp Random.hpp Link.cpp Link.hpp Potentials.cpp Potentials.hpp System.cpp System.hpp worm_canonical.cpp input.example Makefile
