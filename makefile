# Name of the project
PROJ_NAME := MoleKing_util$(shell python3-config --extension-suffix)

OS := $(shell uname)

ifeq ($(OS), Windows_NT)
	CC_FLAGS= -O3 -shared -std=c++11 -I <path-to-pybind11>/include `python-config --cflags --ldflags`
else
	ifeq ($(OS), Darwin)
		CC_FLAGS=-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes`
	else
		CC_FLAGS=-O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes`
	endif
endif
 
all:	
	mkdir MoleKing_util
	mv $(PROJ_NAME) MoleKing_util/
	rm -rf ./src/*/*.o
	rm -rf ./src/*.o
	@echo ''
	@echo "Sucess!"
	@echo "MoleKing_util was compiled for a "$(OS)" System."
		

# .c files
C_SOURCE=$(wildcard  ./src/*/*.cpp)
 
# .h files
H_SOURCE=$(wildcard ./src/*/*.hh)
 
# Object files
OBJ=$(C_SOURCE:.cpp=.o)
 
# Compiler
CC := c++
 
# Compilation and linking
#
all: $(PROJ_NAME)

$(PROJ_NAME): $(OBJ) main.o
	$(CC) -o $@ $^
 
%.o: %.cpp %.hh
	$(CC) -o $@ $< $(CC_FLAGS)
 
main.o: ./src/main.cpp $(H_SOURCE)
	$(CC) $(CC_FLAGS) $<  -o $@ 
 
clean:
	rm -rf ./src/*/*.o
	rm -rf ./src/*.o
	rm -r ./MoleKing_util

