# Name of the project
PROJ_NAME := MoleKing_util$(shell python3-config --extension-suffix)

OS := $(shell uname)

ifeq ($(OS), Windows_NT)
	CC_FLAGS= -O3 -shared -std=c++11 -I <path-to-pybind11>/include `python-config --cflags --ldflags`
else ifeq ($(OS), Darwin)
	CC_FLAGS=-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes`
else
	CC_FLAGS=-O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes`
endif


# .c files
C_SOURCE=$(wildcard  ./src/*/*.cpp)
 
# .h files
H_SOURCE=$(wildcard ./src/*/*.hh)
 
O_PATH= ./bin/

# Object files
OBJ= $(O_PATH)Hessian.o $(O_PATH)AtomicScale.o $(O_PATH)Molecule.o $(O_PATH)PeriodicTable.o $(O_PATH)Geometry.o $(O_PATH)MassCenter.o $(O_PATH)Matrix.o $(O_PATH)Matrix.o

# Compiler
CC := c++
 
# Compilation and linking
#
all: $(PROJ_NAME)
	@echo ''
	@echo "####################################################"
	@echo ''
	@echo "                   Sucess!"
	@echo " MoleKing_util was compiled for a "$(OS)" System."
	@echo ''
	@echo "####################################################"
	@echo ''

$(PROJ_NAME): $(OBJ) $(O_PATH)main.o
	$(CC) -o $@ $^

$(O_PATH)Hessian.o: ./src/berny/Hessian.cpp ./src/berny/Hessian.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)AtomicScale.o: ./src/chemicalUnits/AtomicScale.cpp ./src/chemicalUnits/AtomicScale.hh
	$(CC) -o $@ $< $(CC_FLAGS)
 
$(O_PATH)Molecule.o: ./src/chemicalUnits/Molecule.cpp ./src/chemicalUnits/Molecule.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)PeriodicTable.o: ./src/chemicalUnits/PeriodicTable.cpp ./src/chemicalUnits/PeriodicTable.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)SupraMolecule.o: ./src/chemicalUnits/SupraMolecule.cpp ./src/chemicalUnits/SupraMolecule.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)Geometry.o: ./src/math/Geometry.cpp ./src/math/Geometry.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)MassCenter.o: ./src/math/MassCenter.cpp ./src/math/MassCenter.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)Matrix.o: ./src/math/Matrix.cpp ./src/math/Matrix.hh
	$(CC) -o $@ $< $(CC_FLAGS)

$(O_PATH)main.o: ./src/main.cpp $(H_SOURCE)
	$(CC) $(CC_FLAGS) $<  -o $@ 
 
clean:
	rm -rf bin
	mkdir bin
	mkdir ./bin/include


