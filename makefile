#######  MoleKing_util MakeFile  #######


# Type here
PYTHON3_PATH:= MoleKing_util$(shell python3-config --extension-suffix)


# Compilation and linking
#
all: install
	@echo oi
 
install:
	python3 setup.py

clean:
	rm -rf bin
	mkdir bin
	mkdir ./bin/include


