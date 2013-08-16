run: nanotubo
	./nanotubo

nanotubo:
	gfortran -O3 nanotubo.f90 -o nanotubo

data: nanotubo
	./nanotubo > data
