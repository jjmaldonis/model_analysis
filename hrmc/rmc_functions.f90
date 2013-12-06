!*********************************************************
!**********This module is written by FY on 12/19/2008****
!**********Include three functions***********************
!*********************************************************

!    function chi_square, check_cutoffs and generate_move are included
!

! Changelog
!    unused arguments (r_e etc) removed from chi_sqaure, 12/20/08 pmv
!    subroutine random_move updated. 01/12/08 Jinwoo Hwang
!    function check_cufoff is completed by Feng Yi on 01/16/2009, may need a test program
!    check_cutoffs is modified by Feng Yi on 01/21/2009
!    removed "use ReadInputs" which doesn't seem to be needed; pmv 1/23/09
!    added interface block for function ran2 to enable clean compilation pmv 03/12/09
!    In check_cutoffs, pbc is added to the pair calculation after the hutch is calling. jwh 04/14/2009
!    Added deallocate(atoms) in check_cutoffs, pmv 04/17/09
!    change if condition when deallocate(atoms) in check_cutoffs, FY 04/17/2009
!    gr_e_err is added - jwh 04/25/2009

MODULE rmc_functions

  use model_mod
  implicit none

  integer, dimension(:), allocatable :: total_exp, indx3p_exp, indx4p_exp, indx5p_exp, indx6p_exp
  real :: acceptance_crit
  integer :: vp_we_care_about

  interface
     function ran2(idum)
       real :: ran2
       integer :: idum
     end function ran2
  end interface

CONTAINS

    subroutine setup_vps(natoms)
    integer, intent(in) :: natoms
    integer :: vp_count, istat
    integer :: total, i3, i4, i5, i6, j, i
    integer :: temp, trash1, trash2, trash3
    
    acceptance_crit = 0.01*natoms

    ! This is the file containing the "experimental" vp from the model we want
    ! to replicate. It must have a -1 as the last line
    open(unit=1,file="VP_stat.in",iostat=istat,status='old')
    read(1,*) ! Comment
    temp = 0
    vp_count = 0
    do while(temp .ne. -1)
        read(1,*) temp
        vp_count = vp_count + 1
    enddo
    vp_count = vp_count - 1

    ! Count the number of vp we care about
    vp_we_care_about = 0
    rewind(1)
    read(1,*) ! Comment
    do i=1,vp_count
        read(1,*) j, total, trash1, trash2, trash3, i3, i4, i5, i6
        if(total .ge. acceptance_crit) then
            vp_we_care_about = vp_we_care_about + 1
        endif
    enddo

    allocate(total_exp(vp_we_care_about))
    allocate(indx3p_exp(vp_we_care_about))
    allocate(indx4p_exp(vp_we_care_about))
    allocate(indx5p_exp(vp_we_care_about))
    allocate(indx6p_exp(vp_we_care_about))

    ! Read in the vp we care about
    rewind(1)
    read(1,*) ! Comment
    temp = 0
    do i=1,vp_count
        read(1,*) j, total, trash1, trash2, trash3, i3, i4, i5, i6
        if(total .ge. acceptance_crit) then
            temp = temp + 1
            total_exp(temp) = total
            indx3p_exp(temp) = i3
            indx4p_exp(temp) = i4
            indx5p_exp(temp) = i5
            indx6p_exp(temp) = i6
        endif
    enddo

    write(*,*) "Writing VP we read in:"
    do i=1, vp_we_care_about
        write(*,*)  total_exp(i), indx3p_exp(i), indx4p_exp(i), indx5p_exp(i), indx6p_exp(i)
    enddo
    close(1)
    end subroutine

    FUNCTION chi_square(natoms, ncount, indx3p, indx4p, indx5p, indx6p, iacu)
        integer, intent(in) :: natoms, iacu
        !integer, dimension(:), intent(in) :: ncount, indx3p, indx4p, indx5p, indx6p
        integer, intent(in) :: ncount(:), indx3p(:), indx4p(:), indx5p(:), indx6p(:)
        integer :: j, i
        real :: chi_square, vp_ratio
        chi_square = 0.0

        ! ncount(i) is total of vp_i. indx3p, indx4p, indx5p, indx6p are the vp
        ! indicies. iacu is the end of the list.

        do i=1,iacu
            do j=1, vp_we_care_about
                if( indx3p(i) .eq. indx3p_exp(j) .and. indx4p(i) .eq. indx4p_exp(j) .and. indx5p(i) .eq. indx5p_exp(j) .and. indx6p(i) .eq. indx6p_exp(j) ) then
                    chi_square = chi_square + (total_exp(j) - ncount(i))**2
                    !! We have found the matching vp in the current model.
                    !! See how many we have, and compare that to how many we
                    !! want. Make this a ratio, which we want to be always < 1.
                    !vp_ratio = real(ncount(i))/real(total_exp(j))
                    !if(vp_ratio .gt. 1.0) vp_ratio = 1.0/vp_ratio
                    !! Weight the sum by the amount of vps in the "experiment".
                    !chi_square = chi_square + (vp_ratio * (real(total_exp(j))/real(natoms)))
                endif
            enddo
        enddo

        !! chi_square should be 0 < chi_square < 1.
        !! Higher values represents a better fit.
        !write(*,*) "chi2_vor =", chi_square
    end function chi_square


    function check_cutoffs(m,cutoff_r,moved_atom)
        logical check_cutoffs
        real,dimension(:,:) ::cutoff_r
        integer moved_atom
        type(model) m
        integer, dimension(:), pointer :: atoms
        integer  nlist
        integer istat
        real radius, temp_x, temp_y, temp_z
        integer i,j
        integer num1, num2
        real dist_pair

        !find the maximum cut-off distance
        radius=maxval(maxval(cutoff_r, 1),1)

        ! Find all atoms within radius radius of moved_atom and put them 
        ! in the list atoms. Also sets nlist == size(atoms)+1.
        call hutch_list_3D(m, m%xx%ind(moved_atom),m%yy%ind(moved_atom),m%zz%ind(moved_atom), radius, atoms, istat, nlist)
        if (istat .eq. 1) then
            print *, 'memory allocation fails!'
            return
        else if(istat .eq. -1) then
            print *, 'no atom is found'
            check_cutoffs = .true. 
            return
        endif

        !begin to calculate pair distance
        !note, there are (nlist-1) atoms in total in 'atoms'

        !first, determine the type of moved_atom 
        do i=1, m%nelements
            if(m%znum%ind(moved_atom) .eq. m%atom_type(i)) then
                num1 = i
                exit
            endif
        enddo
        
        do i=1, (nlist-1)
            !do not count the pair distance to itself
            ! The below if statement should be irrelevant. moved_atom isn't put
            ! in atoms by hutch_list_3d as far as i know.
            if (atoms(i) .ne. moved_atom) then
                do j=1, m%nelements
                    if(m%atom_type(j) .eq. m%znum%ind(atoms(i))) then
                        num2 = j
                        exit
                    endif
                enddo !j=1, m%nelements
              
                !endif  !032409 - jwh
               
                !calculate the atomic distance
                !compare with cutoff_r
                temp_x = abs(m%xx%ind(moved_atom) - m%xx%ind(atoms(i)))  !pbc added - jwh 04/14/2009
                temp_y = abs(m%yy%ind(moved_atom) - m%yy%ind(atoms(i)))
                temp_z = abs(m%zz%ind(moved_atom) - m%zz%ind(atoms(i)))

                temp_x = temp_x - m%lx*anint(temp_x/m%lx)
                temp_y = temp_y - m%ly*anint(temp_y/m%ly)
                temp_z = temp_z - m%lz*anint(temp_z/m%lz)
                  
                dist_pair = temp_x**2 + temp_y**2 + temp_z**2
                dist_pair = sqrt(dist_pair)
                   
                if (dist_pair  .lt. cutoff_r(num1, num2)) then
                    check_cutoffs=.false.
                    exit
                else
                    check_cutoffs = .true.  !jwh 032409
                endif
            endif !032409 - jwh
        enddo !i=1, (nlist-1)

        !if (associated(atoms)) deallocate(atoms) ! pmv 4/17/09
        if (nlist .gt. 1) deallocate(atoms) ! fy 4/17/09
    end function check_cutoffs
  


  subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
    use mpi  
    type(model), intent(inout) :: m
    integer iseed, w
    real alpha, aa, bb, cc
    real, intent(out) :: xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new    !Cur and new positions of atoms
    real :: rand1, rand2, rand3, rand4

    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    rand3 = ran2(iseed)
    rand4 = ran2(iseed)

    w = int(m%natoms*rand1)+1

    !write(*,*)myid, rand1, rand2, rand3, rand4

    xx_cur = m%xx%ind(w)            !Cur positions of the atom before random move
    yy_cur = m%yy%ind(w)
    zz_cur = m%zz%ind(w)
    
    !aa = alpha*(rand2 - 0.5)
    !bb = alpha*(rand3 - 0.5)
    !cc = alpha*(rand4 - 0.5)
    aa = alpha*(2*rand2 - 1.0)
    bb = alpha*(2*rand3 - 1.0)
    cc = alpha*(2*rand4 - 1.0)
    
    m%xx%ind(w) = m%xx%ind(w) + aa 
    m%yy%ind(w) = m%yy%ind(w) + bb 
    m%zz%ind(w) = m%zz%ind(w) + cc
    
    if(m%xx%ind(w)>m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)-m%lx       !pbc 
    if(m%yy%ind(w)>m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)-m%ly
    if(m%zz%ind(w)>m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)-m%lz
    if(m%xx%ind(w)<-m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)+m%lx
    if(m%yy%ind(w)<-m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)+m%ly
    if(m%zz%ind(w)<-m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)+m%lz
    
    xx_new=m%xx%ind(w)              !new positions of the atom after random move
    yy_new=m%yy%ind(w)
    zz_new=m%zz%ind(w)
  end subroutine random_move


END MODULE rmc_functions
