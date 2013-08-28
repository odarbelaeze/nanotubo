program nanotubo

    implicit none

    ! Parameters

    double precision, parameter :: pi = 3.14159265359

    ! Auxiliar counters and accumulators

    integer :: i, j, k, err

    ! Data structures
    double precision, allocatable, dimension(:,:) :: r, s
    integer, allocatable, dimension(:) :: nnb
    integer, allocatable, dimension(:,:) :: nbh
    integer :: n

    ! Nanotube properties
    integer :: height = 50
    integer :: nring = 12

    ! Physical parameters
    double precision, dimension(3) :: e_h, e_k
    double precision :: n_h, n_k, j_ex


    ! Setting physical parameters
    e_h = (/ 0.0d0, 0.0d0, 1.0d0 /); n_h = 0.0
    e_k = (/ 0.0d0, 0.0d0, 1.0d0 /); n_k = 0.0
    j_ex = 2.0d0


    ! Program

    call nanotube_sample
    call find_nbh(1.2d0)

    call random_spin

    do i = 1, n, 1
        print*, s(:,i)
    end do

    call clear

contains

    subroutine nanotube_sample
        double precision, parameter :: arc_length = 1.73205080756888
        double precision :: x, y, z, theta, phi, radius

        n = nring * height
        theta = 2.0d0 * pi / nring
        phi = theta / 2.0d0
        radius = arc_length / theta
        z = 0.0d0

        allocate(r(3, n), stat=err)
        if (err /= 0) print *, "r: Allocation request denied"

        do i = 1, height

            phi = 0.0d0
            if ( mod(i / 2, 2) /= 0 ) phi = theta / 2.0d0
            
            do j = 0, nring - 1
                x = radius * cos(j * theta + phi)
                y = radius * sin(j * theta + phi)
                r(:, (i - 1) * nring + j + 1) = (/ x, y, z /)
            end do

            if ( mod(i / 2, 2) == mod((i + 1) / 2, 2) ) then
                z = z + 1.0d0
            else
                z = z + 0.5d0
            end if

        end do


    end subroutine nanotube_sample


    subroutine find_nbh(radius)
        double precision, intent(in) :: radius

        allocate(nnb(n), stat=err)
        if (err /= 0) print *, "nnb: Allocation request denied"

        if ( .not. allocated(r) ) then
            print *, "You need to generate the nanotube first"
            stop "Somehow semantic error."
        end if

        nnb = 0

        do i = 1, n, 1
            do j = 1, n, 1
                if ( sqrt(sum( (r(:, i) - r(:, j)) ** 2 ) ) < radius .and. i /= j ) then
                    nnb(i) = nnb(i) + 1
                end if
            end do
        end do

        allocate(nbh(maxval(nnb), n), stat=err)
        if (err /= 0) print *, "nbh: Allocation request denied"

        nnb = 0

        do i = 1, n, 1
            do j = 1, n, 1
                if ( sqrt(sum( (r(:, i) - r(:, j)) ** 2 ) ) < radius .and. i /= j ) then
                    nnb(i) = nnb(i) + 1
                    nbh(nnb(i), i) = j
                end if
            end do
        end do

    end subroutine find_nbh

    subroutine random_spin

        allocate(s(3,n), stat=err)
        if (err /= 0) print *, "s: Allocation request denied"

        call random_number(s)
        s = s - 0.5
        do i = 1, n, 1
            s(:,i) = s(:,i) / sqrt(sum(s(:,i) ** 2))
        end do

    end subroutine random_spin

    function disturbed_spin(radius)
        real, intent(in) :: radius
        double precision, dimension(3, n) :: disturbed_spin
        
        call random_number(disturbed_spin)
        disturbed_spin = disturbed_spin - 0.5
        disturbed_spin = s + 2.0 * radius * disturbed_spin
        do i = 1, n, 1
            disturbed_spin(:,i) = disturbed_spin(:,i) / sqrt(sum(disturbed_spin(:,i) ** 2))
        end do

    end function disturbed_spin

    function energy(i)
        integer, intent(in) :: i
        double precision :: energy

        energy = n_h  * dot_product(s(:,i), e_h) + &
                 n_k  * dot_product(s(:,i), e_k) ** 2 + &
                 j_ex * dot_product(s(:,i), sum(s(:, nbh(1:nnb(i), i)), dim=2))

    end function energy

    subroutine clear

        if (allocated(r)) deallocate(r, stat=err)
        if (err /= 0) print *, "r: Deallocation request denied"

        if (allocated(nnb)) deallocate(nnb, stat=err)
        if (err /= 0) print *, "nnb: Deallocation request denied"
        
        if (allocated(nbh)) deallocate(nbh, stat=err)
        if (err /= 0) print *, "nbh: Deallocation request denied"

        if (allocated(s)) deallocate(s, stat=err)
        if (err /= 0) print *, "s: Deallocation request denied"

    end subroutine clear

end program nanotubo
