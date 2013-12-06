! Reverse Monte Carlo structural simulation USING VORONOI POLYHEDRA STATS

program rmc

    use omp_lib
    use rmc_global
    use readinputs
    use model_mod
    use rmc_functions
    use eam_mod
    use vor_hrmc
    use shared_data
    implicit none
    include 'mpif.h'
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename, vor_param_fn
    character (len=256) :: outbase
    character (len=256) :: jobID, c, step_str
    character (len=512) :: comment
    character (len=256) :: time_elapsed, output_model_fn, energy_fn, final_model_fn, chi_squared_file, acceptance_rate_fn
    real :: temperature
    real :: max_move
    real, pointer, dimension(:,:) :: cutoff_r 
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: chi2_old, chi2_new, del_chi, beta, chi2_no_energy
    real :: alpha
    real :: sigma_E(20) = 0.0, sigma_chi2(20) = 0.0, sig_E, sig_chi2
    integer :: num_accepted = 0
    integer :: i, j
    integer :: w
    integer :: istat, status2, length
    integer :: rejects
    integer :: step_start
    integer :: iseed2
    !integer :: natoms
    real :: randnum
    real :: te1, te2
    logical :: accepted, use_rmc
    integer :: ipvd, nthr
    doubleprecision :: t0, t1, t2 !timers
    integer, dimension(100) :: acceptance_array
    real :: avg_acceptance

    !------------------- Program setup. -----------------!

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)
    
    nthr = omp_get_max_threads()
    if(myid.eq.0)then
        write(*,*)
        write(*,*) "This is the dev version of rmc!"
        write(*,*)
        write(*,*) "Using", numprocs, "processors."
        write(*,*) "OMP found a max number of threads of", nthr
        write(*,*)
    endif
    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        jobID = '_temp'
    endif

    ! Set input filenames.
    param_filename = 'param_file.in'

    ! Set output filenames.
    outbase = ""
    write(time_elapsed, "(A12)") "time_elapsed"
    time_elapsed = trim(trim(time_elapsed)//jobID)//".txt"
    write(output_model_fn, "(A12)") "model_update"
    output_model_fn = trim(trim(output_model_fn)//jobID)//".txt"
    write(final_model_fn, "(A11)") "model_final"
    final_model_fn = trim(trim(final_model_fn)//jobID)//".txt"
    write(energy_fn, "(A6)") "energy"
    energy_fn = trim(trim(energy_fn)//jobID)//".txt"
    write(chi_squared_file, "(A11)") "chi_squared"
    chi_squared_file = trim(trim(chi_squared_file)//jobID)//".txt"
    write(acceptance_rate_fn, "(A15)") "acceptance_rate"
    acceptance_rate_fn = trim(trim(acceptance_rate_fn)//jobID)//".txt"

    !------------------- Read inputs and initialize. -----------------!

    ! Start timer.
    t0 = omp_get_wtime()

    ! Read input model
    call read_model(param_filename, comment, m, natoms, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
    call read_inputs(param_filename, vor_param_fn, model_filename, temperature, max_move, cutoff_r, step_start, status2)

    call setup_vps(natoms)

    if(myid .eq. 0) write(*,*) "Model filename:", trim(model_filename)

    beta=1./((8.6171e-05)*temperature)
    iseed2 = 104756
    use_rmc = .TRUE.

    call read_eam(m)
    call eam_initial(m,te1); if(myid.eq.0) write(*,*) "Energy = ", te1
    call vor_init(vor_param_fn, model_filename)
    call vor_update()

    !------------------- Start RMC. -----------------!

    t1 = omp_get_wtime()
    open(20, file=param_filename,iostat=istat, status='old')
        rewind(20)
        read(20, '(a80)') comment !read comment from paramfile
        read(20, '(a80)') comment !model filename
        read(20, '(a80)') comment !vor_params filename
        read(20, *) i !step to start at (start at 1 not 0)
        read(20, *) temperature !starting temp with the i above
    close(20)

    if(use_rmc) then ! End here if we only want femsim. Set the variable above.

        alpha = 1.0

        ! Calculate initial chi2
        chi2_no_energy = chi_square(natoms, ncount, indx3p, indx4p, indx5p, indx6p, iacu)
        !chi2_old = chi2_no_energy + te1/te1
        chi2_old = chi2_no_energy + alpha*te1/natoms
        ! Use a function to convert chi2 to a version MC will use correctly
        !chi2_old = 1.0/log(chi2_old/2.0 + 1.0) - 1.0/log(2.0)
        e2 = e1

        t0 = omp_get_wtime()
        if(myid.eq.0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step =      ", i
            write(*,*) "   Chi2_vor =   ", chi2_no_energy
            write(*,*) "   Energy =     ", te1
            write(*,*) "   Temperature =", temperature
            write(*,*)
        endif

        ! Reset time_elapsed, energy_function, chi_squared_file
        if(myid .eq. 0) then
        open(35,file=trim(time_elapsed),form='formatted',status='unknown')
            t1 = omp_get_wtime()
            write(35,*) numprocs, "processors are being used."
            write(35,*) "Step || Time elapsed || Avg time per step || This step's time || Approx time remaining"
        close(35)
        open(34,file=trim(energy_fn),form='formatted',status='unknown')
            write(34,*) "step, energy"
            write(34,*) i, te1
        close(34)
        open(36,file=trim(chi_squared_file),form='formatted',status='unknown')
            write(36,*) "step, chi2, energy, chi2-energy, chi2/energy, alpha weight, beta weight"
            write(36,*) i, chi2_no_energy, te1, chi2_old, abs(chi2_no_energy/te1)
        close(36)
        open(37,file=trim(acceptance_rate_fn),form='formatted',status='replace',access='sequential')
            write(37,*) "step, acceptance rate averaged over last 100 steps"
        close(37)
        endif

        rejects = 0

        ! RMC loop begins.
        do while (i > 0)
            if(myid .eq. 0) write(*,*) "Starting step", i
            t2 = omp_get_wtime()

            call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            ! check_curoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m,cutoff_r,w) )
                ! Check_cutoffs returned false so reset positions and try again.
                m%xx%ind(w) = xx_cur
                m%yy%ind(w) = yy_cur
                m%zz%ind(w) = zz_cur
                call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            end do

            ! Update hutches, data for chi2, and chi2/del_chi
            call hutch_move_atom(m,w,xx_new, yy_new, zz_new)
            call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
            call vor_update_pos(w, xx_new, yy_new, zz_new)
            call vor_update()
            
            chi2_no_energy = chi_square(natoms, ncount, indx3p, indx4p, indx5p, indx6p, iacu)
            !chi2_new = chi2_no_energy + te2/te1 ! Scale with te1
            chi2_new = chi2_no_energy + alpha*te2/natoms
            ! Use a function to convert chi2 to a version MC will use correctly
            !chi2_new = 1.0/log(chi2_new/2.0 + 1.0) - 1.0/log(2.0)
            del_chi = chi2_new - chi2_old
            call mpi_bcast(del_chi, 1, mpi_real, 0, mpi_comm_world, mpierr)
!write(*,*) "chi2_vor, alpha*energy, chi2:", chi2_no_energy, alpha*te2/natoms, chi2_new
!write(*,*) "del_chi = ", del_chi
            randnum = ran2(iseed2)
            ! Test if the move should be accepted or rejected based on del_chi
            if(del_chi <0.0)then
                ! Accept the move
                !num_accepted = num_accepted + 1
                !sigma_E(num_accepted) = alpha*te2/natoms
                !sigma_chi2(num_accepted) = chi2_no_energy
                e1 = e2
                chi2_old = chi2_new
                accepted = .true.
                rejects = 0
                if(myid .eq. 0) write(*,*) "MC move accepted outright."
                if(myid .eq. 0) write(*,*) "Chi2_vor, energy:", chi2_no_energy, te2/natoms
            else
                ! Based on the random number above, even if del_chi is negative, decide
                ! whether to move or not (statistically).
                if(log(1.-randnum) < -del_chi*beta)then
                    ! Accept move
                    num_accepted = num_accepted + 1
                    !sigma_E(num_accepted) = alpha*te2/natoms
                    !sigma_chi2(num_accepted) = chi2_no_energy
                    e1 = e2
                    chi2_old = chi2_new
                    accepted = .true.
                    rejects = 0
                    if(myid .eq. 0) write(*,*) "MC move accepted due to probability." ! del_chi*beta = ", del_chi*beta
                    if(myid .eq. 0) write(*,*) "Chi2_vor, energy:", chi2_no_energy, te2/natoms
                else
                    ! Reject move
                    e2 = e1
                    call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                    call vor_update_pos(w, xx_cur, yy_cur, zz_cur)
                    call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                    accepted = .false.
                    rejects = rejects + 1
                    if(myid .eq. 0) write(*,*) "MC move rejected.", chi2_old, chi2_new
                endif
            endif
            
            write(*,*) "DEL_CHI =", del_chi

            if(rejects .ge. 500) then
                rejects = 0
                temperature = temperature * 2.0
            endif

            !if(num_accepted .eq. 20) then
            !    num_accepted = 0
            !    ! Calculate std devs of both
            !    sig_E = 0.0
            !    sig_chi2 = 0.0
            !    do j=1,size(sigma_E)
            !        sig_E = sig_E + (sigma_E(j) - sum(sigma_E)/size(sigma_E))**2
            !        sig_chi2 = sig_chi2 + (sigma_chi2(j) - sum(sigma_chi2)/size(sigma_chi2))**2
            !    enddo
            !    sig_E = sqrt(sig_E/size(sigma_E))
            !    sig_chi2 = sqrt(sig_chi2/size(sigma_chi2))
            !    !write(*,*) sig_E
            !    !write(*,*) 1.0/real(size(sigma_E))**3, size(sigma_E)*sum(sigma_E**2), sum(sigma_E)**2
            !    !write(*,*) sig_chi2
            !    write(*,*) "OLD alpha, chi2_old:", alpha, chi2_old
            !    alpha = sig_chi2/sig_E
            !    ! We KNOW the current move was accepted, and so we can set
            !    ! chi2_old using the new alpha value here.
            !    chi2_old = chi2_no_energy + alpha*te2/natoms
            !    !if(sig_E .eq. 0.0) sig_E = 1.0E-7
            !    write(*,*) "NEW alpha, chi2_old:", alpha, chi2_old
            !    sigma_E = 0.0
            !    sigma_chi2 = 0.0
            !endif
            
            if(accepted) then
                acceptance_array(mod(i,100)+1) = 1
            else
                acceptance_array(mod(i,100)+1) = 0
            endif
            if(i .gt. 100) avg_acceptance = sum(acceptance_array)/100.0

            ! ---------- Periodically save data. ------------- !
            ! --- The only thing below is saving the data. --- !
            ! ------------------------------------------------ !
            if(myid.eq.0)then
            if(mod(i,1000)==0)then
                ! Write to vor_update

                ! Write to model_update
                write(output_model_fn, "(A12)") "model_update"
                write(step_str,*) i
                output_model_fn = trim(trim(trim(trim(output_model_fn)//jobID)//"_")//step_str)//".xyz"
                open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                    write(33,*)"updated model"
                    write(33,*)m%lx,m%ly,m%lz
                    do j=1,m%natoms
                        write(33,*)m%znum%ind(j), m%xx%ind(j), m%yy%ind(j), m%zz%ind(j)
                    enddo
                    write(33,*)"-1"
                close(33)
            endif
            if(mod(i,1)==0)then
                if(accepted) then
                    ! Write to energy_function
                    open(34,file=trim(energy_fn),form='formatted',status='unknown',access='append')
                        !write(*,*) i, te2
                        write(34,*) i, te2
                    close(34)
                    ! Write chi2 info
                    open(36,file=trim(chi_squared_file),form='formatted',status='unknown',access='append')
                        write(36,*) i, chi2_no_energy, te2, chi2_old, abs(chi2_no_energy/te2)
                    close(36)
                endif
            endif
            if(mod(i,1)==0)then
                ! Write to time_elapsed
                open(35,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
                    t1 = omp_get_wtime()
                    !write(*,*) "Step, time elapsed, temp:", i, t1-t0, temperature
                    write (35,*) i, t1-t0, (t1-t0)/i, t1-t2, (t1-t0)/i * (150000 * log(30/temperature)/log(sqrt(0.7)) - i)
                close(35)
            endif
            if(mod(i,100)==0 .and. i .gt. 100)then
                ! Write to acceptance rate
                open(40,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                    write(40,*) i, avg_acceptance
                close(40)
            endif
            endif ! myid == 0

            ! Every 50,000 steps lower the temp, max_move, and reset beta.
            if(mod(i,10000)==0)then
                temperature = temperature * sqrt(0.7)
                if(myid.eq.0)then
                    write(*,*) "Lowering temp to", temperature, "at step", i
                endif
                max_move = max_move * sqrt(0.94)
                beta=1./((8.6171e-05)*temperature)
            endif
            
            !if(myid .eq. 0) write(*,*) "Finished step", i
            i=i+1
        enddo !RMC do loop

        write(*,*) "Monte Carlo Finished!"

        ! The rmc loop finished. Write final data.
        if(myid.eq.0)then
            ! Write final vor stats

            ! Write final model
            open(unit=55,file=trim(final_model_fn),form='formatted',status='unknown')
            write(55,*)"final model"
            write(55,*)m%lx,m%ly,m%lz
            do i=1,m%natoms
                write(55,*)m%znum%ind(i), m%xx%ind(i), m%yy%ind(i), m%zz%ind(i)
            enddo
            write(55,*)"-1"; close(55)
            ! Write final energy.
            open(56,file=trim(energy_fn),form='formatted', status='unknown',access='append')
            write(56,*) i, te2
            close(56)
            ! Write final time spent.
            open(57,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
            t1 = omp_get_wtime()
            write (57,*) i, t1-t0
            write(57,*) "Finshed.", numprocs, "processors."
            close(57)
        endif
    endif ! Use RMC
    call mpi_finalize(mpierr)

end program rmc
