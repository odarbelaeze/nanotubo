program nanotubo
implicit none

double precision, dimension(:,:), allocatable :: r, s

integer :: n, nring, height
double precision :: rad, sep
double precision :: theta, phi, x, y, z
double precision, parameter :: pi = 3.141592653589793

integer :: i, j, k


nring = 10; height = 50; theta = 2.0 * pi / nring;
n = nring * height
rad = 1.0; sep = 1.0

allocate(r(n, 3), s(n, 3))


! Generating the nanotube sample

do i = 0, height - 1
    do j = 0, nring - 1

        if (mod(i, 3) == 0) then
            phi = 0.0;
        else
            phi = theta / 2.0
        endif

        x = rad * cos(j * theta + phi)
        y = rad * sin(j * theta + phi)
        z = i

        r(i * nring + j + 1, :) = (/ x, y, z /)

    enddo
enddo


! Generating random spin initial configuration

call random_number(s)
s = s - 0.5
do i = 1, n
    s(i, :) = s(i, :) / sqrt(sum(s(i, :) ** 2))
enddo


! Writing some stuf to main screen

do i = 1, n
    write(*,*) s(i, :)
enddo

end program nanotubo
