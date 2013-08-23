nanotubo: nanotubo.f90
	gfortran -O3 -o nanotubo nanotubo.f90

run: nanotubo
	./nanotubo

data: nanotubo
	./nanotubo > data
