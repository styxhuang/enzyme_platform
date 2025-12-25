MODE=RELEASE
CC=gcc
ifeq ($(MODE), DEBUG)
CFLAGS= -DDEBUG -Wall -g -Og -fPIC -std=c++11
else
CFLAGS= -O2 -fPIC -std=c++11
endif

INCLUDE_OPENBABEL= -I${CONDA_PREFIX}/include/openbabel3
INCLUDE_PYTHON= -I${CONDA_PREFIX}/include/python3.7m
INCLUDE_EIGEN= -I${CONDA_PREFIX}/include/eigen3
INCLUDE_PYBIND11= -I${CONDA_PREFIX}/lib/python3.7/site-packages/pybind11/include

INCLUDE= ${INCLUDE_EIGEN} ${INCLUDE_OPENBABEL} ${INCLUDE_PYTHON} ${INCLUDE_PYBIND11}

LINK_PATH_LIB= -L${CONDA_PREFIX}/lib
LINK_PYTHON= -lpython3.7m
LINK_OPENBABEL= -lopenbabel
LINKFLAGS= ${LINK_PATH_LIB} ${LINK_OPENBABEL} ${LINK_PYTHON} -lm -lstdc++

source= $(wildcard vina_docking/*.cpp)
cobjects= $(source:%.cpp=%.o)
asmobjs= vina_docking/_vina_engine.cpython-37m-x86_64-linux-gnu.so

all: $(asmobjs)

clean_all:
	rm $(asmobjs) $(cobjects)

clean:
	rm $(cobjects)

$(asmobjs): $(cobjects)
	$(CC) $^  ${LINKFLAGS} -shared -o $@

$(cobjects): %.o : %.cpp
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@


