CXX=pgc++
CXXFLAGS = -O3 -g -std=c++11 -Minfo=accel -fast -O3 -acc -ta=tesla:fastmath -ta:tesla:fma -ta=tesla:cc60 -DGPU
#============================================================
testcomplex: testcomplex.o
	${CXX} ${CXXFLAGS} -o testcomplex testcomplex.o

testcomplex.o: testcomplex.cpp Makefile
	${CXX} ${CXXFLAGS} -c testcomplex.cpp
#============================================================
clean:
	rm -f ${TARGET} testcomplex.o

