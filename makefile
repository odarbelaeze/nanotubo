nanotubo: nanotubo.f90
	gfortran -O3 nanotubo.f90 -o nanotubo

run: nanotubo
	./nanotubo

data: nanotubo
	./nanotubo > data
