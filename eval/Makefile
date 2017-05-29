CXX = g++
VERSION = 0.01
CXXFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
TARGETS = evaluator
OBJ = evaluator.o

all: ${TARGETS}

${TARGETS}: ${OBJ}
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}

clean:
	rm -f ${OBJ} ${TARGETS} *~


