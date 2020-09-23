# Setup directory paths
INC=$(shell pwd)/include
INCDIR=$(shell pwd)/include
SRCDIR=$(shell pwd)/src
BIN=$(shell pwd)/bin
LIB=$(shell pwd)/lib

# specify compiler
CXX = g++

# root flags
ROOT_CFLAGS = $(shell root-config --cflags) -I${INC} -g -std=c++11 -Ofast 
ROOT_LIBS = $(shell root-config --libs)

# all source code to be compiled .cxx .h

CLASSES = CumulantProfileContainer analysisUtil MomentFunctions InputParameterList EventLoop ProtonEfficiency StFemtoTrack StFemtoEvent QALoop MomentEventLoop SystematicEventLoop FlucContainer PileUpCorrection

DIR_H = $(addprefix ${INC}/,${CLASSES})
H_FILES = $(addsuffix .h,${DIR_H})

$(info H_FILES are [${H_FILES}] )

OBJ_FILES = $(addsuffix .o,${CLASSES})

OBJ_FULLPATH = $(addprefix ${BIN}/,${OBJ_FILES})

all: ${LIB}/momentCode.so

${LIB}/momentCode.so: ${BIN}/dict.o ${OBJ_FULLPATH}
	${CXX} -shared -o $@ $^ ${ROOT_LIBS}

${BIN}/dict.C: ${H_FILES}
	rootcint -f $@ -c $^ ${INC}/LinkDef.h

${BIN}/dict.o: ${BIN}/dict.C
	${CXX} ${ROOT_CFLAGS} -fPIC -c $< -o $@

${BIN}/%.o: ${SRCDIR}/%.cxx ${INCDIR}/%.h
	${CXX} ${ROOT_CFLAGS} -fPIC -c $< -o $@


clean:
	rm -rf ${BIN}/*.o ${LIB}/*.so ${BIN}/dict*
