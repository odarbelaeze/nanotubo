program nanotubo
implicit none

double precision, dimension(:,:), allocatable :: r, s

integer :: n, nring, height
double precision :: rad, sep

integer :: i, j, k


nring = 10; height = 50;
rad = 1.0; sep = 1.0

! Generating the nanotube sample

call sample(r, n, height, nring, rad, sep)


! Generating random spin initial configuration

call random_spin(s, n)


! Writing some stuf to main screen

do i = 1, n
    write(*,*) s(i, :)
enddo

contains

subroutine sample(r, n, height, nring, rad, sep)

    double precision, dimension(:,:), allocatable, intent(out) :: r
    integer, intent(out) :: n

    double precision, intent(in) :: rad, sep
    integer, intent(in) :: nring, height

    double precision :: theta, phi, x, y, z
    double precision, parameter :: pi = 3.141592653589793
    integer :: i, j, k

    theta = 2.0 * pi / nring
    n = nring * height

    allocate(r(n, 3))

    do i = 0, height - 1
        do j = 0, nring - 1

            if (mod(i, 3) == 0) then
                phi = 0.0;
            else
                phi = theta / 2.0
            endif

            ! phi = mod(i, 3) == 0? 0.0 : theta / 2.0

            x = rad * cos(j * theta + phi)
            y = rad * sin(j * theta + phi)
            z = i

            r(i * nring + j + 1, :) = (/ x, y, z /)

        enddo
    enddo

end subroutine


subroutine random_spin(s, n)
    double precision, dimension(:,:), allocatable, intent(inout) :: s
    integer, intent(in) :: n

    allocate(s(n, 3))

    call random_number(s)
    s = s - 0.5
    do i = 1, n
        s(i, :) = s(i, :) / sqrt(sum(s(i, :) ** 2))
    enddo

end subroutine

end program nanotubo



