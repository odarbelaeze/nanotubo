program nanotubo

    use utils

    implicit none

    double precision, allocatable, dimension(:,:) :: r, s
    integer, allocatable, dimension(:) :: nnb
    integer, allocatable, dimension(:,:) :: nbh

    integer :: n, nring, height
    double precision :: rad, sep

    integer :: i, j, k


    nring = 10; height = 50;
    rad = 1.0; sep = 1.0

    ! Generating the nanotube sample

    call nanotube_sample(r, n, height, nring, rad, sep)
    call find_nbh(r, 2.25d0, nnb, nbh)


    ! Generating random spin initial configuration

    call random_spin(s, n)

    ! Writing some stuf to main screen

    do i = 1, 10000, 1
        call mc_steep(r, s, nnb, nbh, 2.0d0, 0.0d0, (/ 0.0d0, 0.0d0, 0.0d0 /), 0.01d0)
        write(*,*) i, total_energy(r, s, nnb, nbh, 1.0d0, 0.0d0, (/ 0.0d0, 0.0d0, 0.0d0 /)), sqrt(sum((sum(s, 1) / n) ** 2))
    end do

end program nanotubo



