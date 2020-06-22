CXX=pgc++
CXXFLAGS = -O3 -g -std=c++11 -Minfo=accel -fast -O3 -acc -ta=tesla:fastmath -ta:tesla:fma -ta=tesla:cc60 -DGPU
#============================================================
testreal: testreal.o
	${CXX} ${CXXFLAGS} -o testreal testreal.o

testreal.o: testreal.cpp Makefile
	${CXX} ${CXXFLAGS} -c testreal.cpp
#============================================================
clean:
	rm -f ${TARGET} testreal.o


