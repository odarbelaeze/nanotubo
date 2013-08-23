program nanotubo

    use utils

    implicit none

    double precision, allocatable, dimension(:,:) :: r, s
    integer, allocatable, dimension(:) :: nnb
    integer, allocatable, dimension(:,:) :: nbh

    integer :: n, nring, height
    double precision :: rad, sep

    double precision, allocatable, dimension(:) :: energy_sample, mag_sample
    double precision :: t_step, t_max, temperature
    integer :: n_temperatures, sample_steps, relax_steps

    integer :: i, j, k, err


    nring = 10; height = 51;
    rad = 1.0; sep = 1.0

    ! Generating the nanotube sample

    call nanotube_sample(r, n, height, nring, rad, sep)
    call find_nbh(r, 3.0d0, nnb, nbh)

    print*, "#", n, maxval(nnb, dim=1)

    ! Generating random spin initial configuration

    call random_spin(s, n)

    t_max = 1.5
    t_step = 0.005
    n_temperatures = t_max / t_step

    print*, "#", n_temperatures

    relax_steps = 2000
    sample_steps = 1

    allocate(energy_sample(sample_steps), stat=err)
    if (err /= 0) print *, "# energy_sample: Allocation request denied"

    allocate(mag_sample(sample_steps), stat=err)
    if (err /= 0) print *, "# mag_sample: Allocation request denied"

    do i = 1, n_temperatures
        temperature = t_max - (i - 1) * t_step

        do j = 1, relax_steps
            call mc_steep(r, s, nnb, nbh, 1.0d0, 5.0d0, (/ 0.0d0, 0.0d0, 1.0d0 /), temperature)
        end do

        do j = 1, sample_steps
            call mc_steep(r, s, nnb, nbh, 1.0d0, 5.0d0, (/ 0.0d0, 0.0d0, 1.0d0 /), temperature)
            energy_sample(j) = total_energy(r, s, nnb, nbh, 2.0d0, 1.0d0, (/ 0.0d0, 0.0d0, 1.0d0 /))
            mag_sample(j) = sqrt(sum((sum(s, dim=1) / n) ** 2))
        end do

        print*, temperature, mean(energy_sample, sample_steps), mean(mag_sample, sample_steps), &
                variance(energy_sample, sample_steps)
    
    end do

    if (allocated(energy_sample)) deallocate(energy_sample, stat=err)
    if (err /= 0) print *, "# energy_sample: Deallocation request denied"

    if (allocated(mag_sample)) deallocate(mag_sample, stat=err)
    if (err /= 0) print *, "# mag_sample: Deallocation request denied"

end program nanotubo
