nanotubo: nanotubo.f90  utils.o
	gfortran -O3 -o nanotubo nanotubo.f90 utils.o

utils.o: utils.f90
	gfortran -O3 -c utils.f90

run: nanotubo
	./nanotubo

data: nanotubo
	./nanotubo > data
