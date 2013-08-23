module utils

    implicit none
    contains

    subroutine nanotube_sample(r, n, height, nring, rad, sep)

        double precision, dimension(:,:), allocatable, intent(out) :: r
        integer, intent(out) :: n

        double precision, intent(in) :: rad, sep
        integer, intent(in) :: nring, height

        double precision :: theta, phi, x, y, z
        double precision, parameter :: pi = 3.141592653589793
        integer :: i, j

        theta = 2.0 * pi / nring
        n = nring * height

        allocate(r(n, 3))

        z = 0.0

        do i = 0, height - 1
            do j = 0, nring - 1

                if (mod(i, 3) == 0) then
                    phi = 0.0;
                else
                    phi = theta / 2.0
                endif

                x = rad * cos(j * theta + phi)
                y = rad * sin(j * theta + phi)

                r(i * nring + j + 1, :) = (/ x, y, z /)

                if ( mod(i, 3) /= 0 .and. mod(i + 1, 3) /= 0 ) then
                    z = z + 1.0
                end if

                z = z + 1.0

            enddo
        enddo

    end subroutine


    subroutine find_nbh(r, rad, nnb, nbh)
        double precision, dimension(:,:), allocatable, intent(in) :: r
        double precision, intent(in) :: rad
        integer, allocatable, dimension(:), intent(out) :: nnb
        integer, allocatable, dimension(:,:), intent(out) :: nbh

        integer, dimension(2) :: n
        integer :: i, j, cont

        n = shape(r);
        if ( .not. allocated(nnb) ) allocate(nnb(n(1)))
        nnb = 0

        do i = 1, n(1)
            do j = 1, n(1)
                if ( (sqrt(sum((r(j,:) - r(i,:)) ** 2)) <= rad) .and. (i /= j) ) then
                    nnb(i) = nnb(i) + 1
                end if
            end do
        end do

        if (allocated(nbh)) deallocate(nbh)
        allocate(nbh(n(1), maxval(nnb)))

        do i = 1, n(1)
            cont = 1
            do j = 1, n(1)
                if ( (sqrt(sum((r(j,:) - r(i,:)) ** 2)) <= rad) .and. (i /= j) ) then
                    nbh(i, cont) = j
                    cont = cont + 1
                end if
            end do
        end do

    end subroutine find_nbh


    subroutine random_spin(s, n)
        double precision, dimension(:,:), allocatable, intent(out) :: s
        integer, intent(in) :: n

        integer :: i

        allocate(s(n, 3))

        call random_number(s)
        s = s - 0.5
        do i = 1, n
            s(i, :) = s(i, :) / sqrt(sum(s(i, :) ** 2))
        enddo

    end subroutine


    function energy(i, r, s, nnb, nbh, j_ex, k_an, e_an)
        integer, intent(in) :: i
        double precision, dimension(:,:), allocatable, intent(in) :: r, s
        integer, allocatable, dimension(:), intent(in) :: nnb
        integer, allocatable, dimension(:,:), intent(in) :: nbh

        double precision, intent(in) :: j_ex
        double precision, intent(in) :: k_an
        double precision, dimension(3), intent(in) :: e_an

        double precision :: energy

        integer :: j

        energy = 0

        do j = 1, nnb(i)
            energy = energy - 0.5 * j_ex * dot_product(s(i,:), s(nbh(i,j),:))
        end do

        energy = energy - k_an * dot_product(s(i,:),e_an) ** 2

    end function energy


    function rel_energy(i, r, s, nnb, nbh, j_ex, k_an, e_an)
        integer, intent(in) :: i
        double precision, dimension(:,:), allocatable, intent(in) :: r, s
        integer, allocatable, dimension(:), intent(in) :: nnb
        integer, allocatable, dimension(:,:), intent(in) :: nbh

        double precision, intent(in) :: j_ex
        double precision, intent(in) :: k_an
        double precision, dimension(3), intent(in) :: e_an

        double precision :: rel_energy

        integer :: j

        rel_energy = 0

        do j = 1, nnb(i)
            rel_energy = rel_energy - j_ex * dot_product(s(i,:), s(nbh(i,j),:))
        end do

        rel_energy = rel_energy - k_an * dot_product(s(i,:),e_an) ** 2

    end function rel_energy


    function total_energy(r, s, nnb, nbh, j_ex, k_an, e_an)
        double precision, dimension(:,:), allocatable, intent(in) :: r, s
        integer, allocatable, dimension(:), intent(in) :: nnb
        integer, allocatable, dimension(:,:), intent(in) :: nbh

        double precision, intent(in) :: j_ex
        double precision, intent(in) :: k_an
        double precision, dimension(3), intent(in) :: e_an

        double precision :: total_energy

        integer :: i, j
        integer, dimension(2) :: n
        n = shape(r);

        total_energy = 0

        do i = 1, n(1)
            total_energy = total_energy + energy(i, r, s, nnb, nbh, j_ex, k_an, e_an)
        end do

    end function total_energy

    subroutine mc_steep(r, s, nnb, nbh, j_ex, k_an, e_an, betainv)
        double precision, dimension(:,:), allocatable, intent(in) :: r
        double precision, dimension(:,:), allocatable, intent(inout) :: s
        integer, allocatable, dimension(:), intent(in) :: nnb
        integer, allocatable, dimension(:,:), intent(in) :: nbh
        double precision, intent(in) :: j_ex
        double precision, intent(in) :: k_an
        double precision, dimension(3), intent(in) :: e_an
        double precision, intent(in) :: betainv

        double precision, dimension(:,:), allocatable :: st

        double precision :: e_inicial, e_final, delta_e, random
        double precision, dimension(3) :: s_old

        integer :: i, err = 0
        integer, dimension(2) :: n
        n = shape(r);

        call random_spin(st, n(1));

        do i = 1, n(1)
            e_inicial = rel_energy(i, r, s, nnb, nbh, j_ex, k_an, e_an)
            s_old = s(i, :)
            s(i, :) = st(i, :)
            e_final = rel_energy(i, r, s, nnb, nbh, j_ex, k_an, e_an)
            delta_e = e_final - e_inicial

            if (delta_e > 0) then
                call random_number(random)
                if ( exp( - delta_e / betainv) < random ) then
                    s(i, :) = s_old
                end if
            end if

        end do

        if (allocated(st)) deallocate(st, stat=err)
        if (err /= 0) print *, "# st: Deallocation request denied"

    end subroutine mc_steep

    function mean(array, n)
        integer :: n
        double precision, dimension(n), intent(in) :: array
        double precision :: mean

        mean = sum(array) / n

    end function mean


    function variance(array, n)
        integer :: n
        double precision, dimension(n), intent(in) :: array
        double precision :: variance

        variance = (mean(array ** 2, n)) - (mean(array, n) ** 2)

    end function variance

end module utils
