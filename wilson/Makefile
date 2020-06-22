#============================================================
TARGET  = 2D-Wilson-accel
SOURCES = main.cpp
OBJS    = main.o
INC_PATH=../include
 
CXX=pgc++
CXXFLAGS = -O3 -g -std=c++11 -Minfo=accel -fast -O3 -acc -ta=tesla:fastmath -ta:tesla:fma -ta=tesla:cc60 -DGPU -I${INC_PATH}
#============================================================
 
all: ${TARGET}

${TARGET}: main.o
	$(CXX) ${CXXFLAGS} -o ${TARGET} main.o

main.o: main.cpp Makefile
	${CXX} ${CXXFLAGS} -c main.cpp

#============================================================

clean:
	rm -f ${TARGET} main.o

