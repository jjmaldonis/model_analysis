
module ReadInputs

contains
    subroutine read_inputs(param_filename, vor_param_fn, model_fn, temperature, max_move, cutoff_r, &
        step_start, status2)
        implicit none
        character (len=*), intent(inout) :: param_filename, vor_param_fn !Assume size array, be careful
        character (len=80), intent(out) :: model_fn
        real, intent(out) :: temperature
        real, intent(out) :: max_move
        real, intent(out), dimension(:,:) :: cutoff_r
        integer, intent(out) :: step_start
        integer, intent(out) :: status2

        ! Local variables
        character (len=80) comment1  

        open(20, file=param_filename,iostat=status2, status='old')
        rewind(20)
        if(status2 .ne. 0) then !open fails
            print *, 'cannot open file with name: ', param_filename
            return
        endif
        read(20, '(a80)') comment1 ! read comment and it is neglected in the input
        read(20, '(a80)') model_fn; model_fn = adjustl(model_fn)
        read(20, '(a80)') vor_param_fn; vor_param_fn = adjustl(vor_param_fn)
        read(20, *) step_start! starting step #
        read(20, *) temperature
        read(20, *) max_move
        read(20, * ) !commented line
        read(20, *) cutoff_r ! This is a nelements by nelements matrix

        !********************************************************************
        !*******************READ Voronoi Statistics file*********************
        !********************************************************************
        close(20)

    end subroutine read_inputs

end module readinputs
