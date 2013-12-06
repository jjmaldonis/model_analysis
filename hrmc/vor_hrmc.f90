

module shared_data
    implicit none
    real, parameter :: n0 = 10000
    real, parameter :: maxcan = 75 !75
    real, parameter :: maxver = 100 !100
    real, parameter :: maxepf = 20  !20
    real, parameter :: atol = 0.01
    real, parameter :: tol = 0.01
    integer, save :: nsp, atomtype_num
    integer, save :: ndata, natoms !nsp1, nsp2, nsp3
    integer, save, dimension(:), allocatable :: sp_natoms
    integer, save, dimension(:), allocatable :: znum_sp
    integer, save, dimension(:), allocatable :: id
    real, save :: world_size, rcut
    real, save, dimension(:,:), allocatable :: pos     ! atomic coordinates
    real, save, dimension(:), allocatable :: radii_sp
    real, save :: scx(n0,25), scy(n0,25), scz(n0,25) ! atomic coordinate of the cluster
    real, save :: sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of the center atom at cluster
    real, save :: vvol(n0)  ! vor cell volume
    integer, save :: mscce(n0), msce(n0,25), ibad(n0), mvijk(maxver,3)
    integer, save :: indx3(n0),indx4(n0),indx5(n0),indx6(n0), indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
    integer, save :: nc, nv, ne, nf ! # of candidates,vertices,edges,faces
    real, save :: p(maxcan,4),v(maxver,3) ! vertices
    integer, save :: mtag(maxcan) ! tag of atom after sorting
    integer, save :: nloop(maxepf,maxcan),nepf(maxcan) ! check connection
    real, save :: area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
    integer, save :: nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
    integer, save, dimension(:,:), allocatable ::  nnabsp ! Number of neighbor species???
    integer, save, dimension(:), allocatable :: indx30,indx40,indx50,indx60, indx70,indx80,indx90,indx100,nbr10,nbr20,nbr30
    integer, save, dimension(:), allocatable :: indx3p,indx4p,indx5p,indx6p, indx7p,indx8p,indx9p,indx10p,nbr1p,nbr2p,nbr3p
    integer, save, dimension(:), allocatable :: ncount,ncount1,ncount2,ncount3,ntype
    integer, save :: iacu, iflg, noi, nbad

end module shared_data


module vor_hrmc

contains

function connect(i,j)
    use shared_data
    implicit none
    integer, intent(in) :: i, j
    integer np, ii, jj
    logical :: connect
    np=0
    do ii=1,3
        do jj=1,3
            if(mvijk(i,ii) .eq. mvijk(j,jj)) then
                np = np + 1
            endif
        enddo
    enddo

    select case(np)
    case (0)
        stop 'not connected, same ic?'
    case (1)
        connect=.false.
    case (2)
        connect=.true.
    case (3)
        stop 'wrong connection'
    end select
    return
end function connect


subroutine vor_init(param_filename, model_filename)
    use shared_data
    implicit none
    ! Input parameters
    character(len=256), intent(inout) :: param_filename, model_filename
    character(len=256) :: comment
    integer :: istat
    ! Loop counters
    integer i, k, ii

    ! Open files.
    open(unit=4,file=param_filename, iostat=istat)
    if(istat .ne. 0) then
        write(*,*) "Error opening parameters ilfe, "," ",param_filename," status=",istat
        stop
    endif

    read(4,*)
    read(4,*) nsp
    allocate(sp_natoms(nsp))
    read(4,*)
    read(4,*) ndata, natoms, sp_natoms
    read(4,*)
    read(4,*) world_size, rcut
    close(4)

    allocate(pos(natoms,3))
    allocate(znum_sp(nsp)); znum_sp = 0
    allocate(radii_sp(nsp)); radii_sp = 0
    allocate(id(natoms))
    allocate(nnabsp(natoms,nsp)) ! natoms should be n0 but we get a compile error; natoms should work tho

    allocate(indx30(ndata*natoms))
    allocate(indx40(ndata*natoms))
    allocate(indx50(ndata*natoms))
    allocate(indx60(ndata*natoms))
    allocate(indx70(ndata*natoms))
    allocate(indx80(ndata*natoms))
    allocate(indx90(ndata*natoms))
    allocate(indx100(ndata*natoms))
    allocate(nbr10(ndata*natoms))
    allocate(nbr20(ndata*natoms))
    allocate(nbr30(ndata*natoms))
    allocate(indx3p(ndata*natoms))
    allocate(indx4p(ndata*natoms))
    allocate(indx5p(ndata*natoms))
    allocate(indx6p(ndata*natoms))
    allocate(indx7p(ndata*natoms))
    allocate(indx8p(ndata*natoms))
    allocate(indx9p(ndata*natoms))
    allocate(indx10p(ndata*natoms))
    allocate(nbr1p(ndata*natoms))
    allocate(nbr2p(ndata*natoms))
    allocate(nbr3p(ndata*natoms))
    allocate(ncount(ndata*natoms))
    allocate(ncount1(ndata*natoms))
    allocate(ncount2(ndata*natoms))
    allocate(ncount3(ndata*natoms))
    allocate(ntype(ndata*natoms))

    model_filename = trim(adjustl(model_filename))
    !open(unit=2,file=model_filename,iostat=istat,status='old')
    open(unit=2,file="Zr54Cu38Al8_1523atoms_start.xyz",iostat=istat,status='old')
    rewind(2)
    if(istat .ne. 0) then
        write(*,*) "Error opening model file, "," ",trim(model_filename)," status=",istat
        stop
    endif

    write(*,*) "Assuming", natoms, "atoms in model."
    ! Read in atoms from XYZ file
    read(2,*) comment ! skip comment
    read(2,*) ii ! skip model size line
    atomtype_num = 0
    do ii = 1, natoms
        read(2,*) id(ii), pos(ii,1),pos(ii,2),pos(ii,3)   ! Watch input model format -JWH
        ! Find out if we have a new atom type using this loop
        k = 0
        do i = 1, nsp
            if(id(ii) .eq. znum_sp(i)) k = i
        enddo
        if(k .eq. 0) then ! Then we have a new atom type
            atomtype_num = atomtype_num + 1
            k = atomtype_num
            znum_sp(k) = id(ii) ! Add the new atom to znum_sp; currently id(ii) is the atomic number of the atom
        endif
        id(ii) = k ! Set id(ii) to the right value, cooresponding to the index in znum_sp
    enddo
    close(2)

    ! Fix atomic coordinates to be between 0 and 1.
    if(minval(pos) .ge. -1 .and. maxval(pos) .le. 1 .and. minval(pos) .lt. 0) then ! If fractional coords
        write(*,*) "Started with frac coords"
        pos=0.5*(pos+1.0)
    else if(minval(pos) .lt. -1 .and. maxval(pos) .ge. 1) then ! If normal xyz coords
        write(*,*) "Started with full coords"
        pos = (pos + world_size/2.0)/world_size
    else if(minval(pos) .eq. 0.0 .and. maxval(pos) .gt. 1) then
        pos = pos/maxval(pos)
    else
        write(*,*) "ERROR! Unknown input coordinates. Check input model."
        stop
    endif
end subroutine vor_init


subroutine vor_reset
    use shared_data
    radii_sp = 1.0 ! This is 1 x nsp array

    ! Initiallization; I repalced call initsys with this - Jason
    ! These are all arrays
    nnab = 0
    vvol = 0.0
    nablst = 0
    nedges = 0
    scx = 0.0
    scy = 0.0
    scz = 0.0
    sccx= 0.0
    sccy= 0.0
    sccz= 0.0
    vvol =0.0
    mscce = 0
    msce = 0
    ibad = 0
    mvijk = 0
    indx3  = 0
    indx4  = 0
    indx5  = 0
    indx6  = 0
    indx7  = 0
    indx8  = 0
    indx9  = 0
    indx10  = 0
    nc = 0
    nv = 0
    ne = 0
    nf = 0
    p = 0.0
    v = 0.0
    mtag = 0
    nloop = 0
    nepf = 0
    area = 0.0
    sleng = 0.0
    tleng = 0.0
    nnab = 0
    nablst = 0
    nedges = 0
    nnabsp = 0
    indx30 = 0
    indx40 = 0
    indx50 = 0
    indx60 = 0
    indx70 = 0
    indx80 = 0
    indx90 = 0
    indx100 = 0
    nbr10 = 0
    nbr20 = 0
    nbr30 = 0
    indx3p = 0
    indx4p = 0
    indx5p = 0
    indx6p = 0
    indx7p = 0
    indx8p = 0
    indx9p = 0
    indx10p = 0
    nbr1p = 0
    nbr2p = 0
    nbr3p = 0
    ncount = 0
    ncount1 = 0
    ncount2 = 0
    ncount3 = 0
    ntype = 0
    iacu = 0
    iflg = 0
    noi = 0
    nbad = 0
end subroutine vor_reset


subroutine vor_update_pos(atom, xx, yy, zz)
    use shared_data
    implicit none
    real, intent(in) :: xx, yy, zz
    integer, intent(in) :: atom
    pos(atom,1) = (xx + world_size/2.0)/world_size
    pos(atom,2) = (yy + world_size/2.0)/world_size
    pos(atom,3) = (zz + world_size/2.0)/world_size
end subroutine


subroutine vor_update
    use shared_data
    integer :: i, ii, j, istat
    character(len=256) :: outbase, outfile

    outfile = "vor.out"
    outbase = ""

    call vor_reset

    call vtanal!(maxcan, maxver, maxepf, natoms, world_size, noi, nc, nf, ne, nv, nbad, atol, tol, rcut, id, radii_sp, mtag, nepf, nnab, mscce, ibad, mvijk, nloop, nedges, nablst, msce, area, tleng, vvol, sccx, sccy, sccz, v, pos, p, sleng, scx, scy, scz)
    call outvt!(natoms, nsp, id, nedges, nablst, nnabsp, nnab, indx3, indx4, indx5, indx6, indx7, indx8, indx9, indx10, vvol)

    i=0
    do ii=1, natoms
        i = i + 1
        indx30(ii)=indx3(i)
        indx40(ii)=indx4(i)
        indx50(ii)=indx5(i)
        indx60(ii)=indx6(i)
        indx70(ii)=indx7(i)
        indx80(ii)=indx8(i)
        indx90(ii)=indx9(i)
        indx100(ii)=indx10(i)
        nbr10(ii)=nnabsp(i,1)
        nbr20(ii)=nnabsp(i,2)
        nbr30(ii)=nnabsp(i,3)

        !... ! Figure out how to do this part
        ntype(i) = id(ii)
    enddo

    ! Do the statistics
    iacu=0
    do i=1,ndata*natoms
        iflg=0
        do j=1,iacu
            if(indx30(i).eq.indx3p(j).and.indx40(i).eq.indx4p(j).and.indx50(i).eq.indx5p(j).and.indx60(i).eq.indx6p(j).and.indx70(i).eq.indx7p(j).and. indx80(i).eq.indx8p(j).and.indx90(i).eq.indx9p(j).and.indx100(i).eq.indx10p(j))then
                iflg=1
                ncount(j)=ncount(j)+1
                nbr1p(j)=nbr1p(j)+nbr10(i)
                nbr2p(j)=nbr2p(j)+nbr20(i)
                nbr3p(j)=nbr3p(j)+nbr30(i)

                ! WARNING This only works for ternary or lower I think... TODO
                if(ntype(i).eq.1) then
                    ncount1(j)=ncount1(j)+1
                else
                    if(ntype(i).eq.2) then
                        ncount2(j)=ncount2(j)+1
                    else
                        ncount3(j)=ncount3(j)+1
                    endif
                endif
            endif
        enddo
        if(iflg.eq.0)then
            iacu=iacu+1
            ncount(iacu)=1

            ! WARNING This only works for ternary or lower I think... TODO
            if(ntype(i).eq.1) then
                ncount1(iacu)=1
            else
                if(ntype(i).eq.2) then
                    ncount2(iacu)=1
                else
                    ncount3(iacu)=1
                endif
            endif

            nbr1p(iacu)=nbr10(i)
            nbr2p(iacu)=nbr20(i)
            nbr3p(iacu)=nbr30(i)
            indx3p(iacu)=indx30(i)
            indx4p(iacu)=indx40(i)
            indx5p(iacu)=indx50(i)
            indx6p(iacu)=indx60(i)
            indx7p(iacu)=indx70(i)
            indx8p(iacu)=indx80(i)
            indx9p(iacu)=indx90(i)
            indx10p(iacu)=indx100(i)
        endif
    enddo

    outfile = trim(outbase) // '_stat.out'
    open(unit=8,file=outfile,status='replace', iostat=istat)
    if(istat .ne. 0) then
        write (*,*) "Error opening stat output file, "," ",outfile, " status=",istat
        stop
    endif

    111  format(5I8,8I3,3I10)
    write(8,*) "No., total, center1, center2, center3, n3, n4, n5, n6, n7, n8, n9, n10, nbr-1, nbr-2, nbr-3"
    do i=1,iacu
        write(8,111) i, ncount(i), ncount1(i), ncount2(i),ncount3(i), indx3p(i),indx4p(i),indx5p(i),indx6p(i), indx7p(i),indx8p(i),indx9p(i),indx10p(i), nbr1p(i),nbr2p(i),nbr3p(i)
    enddo

    close(8)
end subroutine vor_update


subroutine vtanal!(maxcan, maxver, maxepf, natoms, world_size, noi, nc, nf, ne, nv, nbad, atol, tol, rcut, id, radii_sp, mtag, nepf, nnab, mscce, ibad, mvijk, nloop, nedges, nablst, msce, area, tleng, vvol, sccx, sccy, sccz, v, pos, p, sleng, scx, scy, scz)
    use shared_data
    implicit none
    !integer, intent(in) :: natoms, maxcan, maxver, maxepf, world_size
    !integer, intent(inout) :: noi, nc, nf, ne, nv, nbad
    !real, intent(in) :: atol, rcut, tol
    !integer, dimension(:), intent(in) :: id, radii_sp
    !integer, dimension(:), intent(inout) :: mtag, nepf, nnab, mscce, ibad
    !integer, dimension(:,:), intent(inout) :: mvijk, nloop, nedges, nablst, msce
    !real, dimension(:), intent(inout) :: area, tleng, vvol, sccx, sccy, sccz
    !real, dimension(:,:), intent(in) :: v
    !real, dimension(:,:), intent(inout) :: pos, p, sleng, scx, scy, scz

    ! Locals
    real :: rxij, ryij, rzij, rijsq, rcutsq, ratio, xij, yij, zij, nx, ny, nz, vol, tarea, sumvol, avgarea, volratio, avglen, x1, x2, y1, y2, z1, z2 !  I think nx, ny, nz are supposed to be reals at least - Jason.
    integer :: i, j, ic, ie, iv, i1, iv1, je, jv, ive, ivs, iicc, iic
    integer, dimension(:) :: nedges1(maxcan)
    !logical :: connect ! Function

    ! For each atom pair, store the distance between them. Tag them too (??).
    sumvol = 0.0
    volratio = 0.0
    vol = world_size**3
    do i=1,natoms
!write(*,*) "i=", i
        ic = 0
        noi = i
        do j=1,natoms
            if(j .ne. i) then
                rxij = pos(j,1)-pos(i,1)
                ryij = pos(j,2)-pos(i,2)
                rzij = pos(j,3)-pos(i,3)

                rxij = rxij - anint ( rxij ) ! JASON I think this is pbc
                ryij = ryij - anint ( ryij )
                rzij = rzij - anint ( rzij )

                ! This line implements weighted voronoi analysis
                ratio=2*radii_sp(id(i))/(radii_sp(id(i))+radii_sp(id(j))) !With all radii = 0.0, the ratio is always zero. -JWH. With all radii = 1.0, the ratio is always 1.0. This doesn't need to be in this loop. - Jason

                rxij=rxij*ratio*world_size
                ryij=ryij*ratio*world_size
                rzij=rzij*ratio*world_size
                rijsq = rxij**2 + ryij**2 + rzij**2
                rcutsq=rcut**2
                ! The below line sets cutoff for each species; optional
                !rcutsq=(rcut*2*rsp(id(i))/(rsp(1)+rsp(2)))**2 ! Commented by Hao and JWH

                if ( rijsq .lt. rcutsq ) then
                    !write(*,*) pos(i,1), pos(i,2), pos(i,3), pos(j,1), pos(j,2), pos(j,3), rijsq
                    ic = ic + 1
                    if(ic.gt.maxcan)then
                        write(*,*)ic, maxcan
                        stop 'too many candidates'
                    endif
                    p(ic,1) = rxij
                    p(ic,2) = ryij
                    p(ic,3) = rzij
                    p(ic,4) = rijsq
                    mtag(ic) = j
                endif
            endif
        enddo

        ! Candidates have been selected (huh??)
        nc = ic
        !write(*,*) "Number of candidates = ", nc
        ! Sort into ascending order of distance
        call sort!(nc, p, mtag)
        call work!(maxver, maxepf, noi, nc, tol, p, v, nepf, nloop, mvijk, nv, nf, ne, nbad, ibad)

        !***************************************************************
        !             sort vertices ring in faces and
        !          calculate edge length and face are
        !***************************************************************

        tarea = 0.0
        avglen = 0.0
        avgarea = 0.0
        ! The below are arrays
        sleng = 0.0
        tleng = 0.0
        area = 0.0

        do ic=1, nc
            if(nepf(ic) .ne. 0) then
                do ie = 1, nepf(ic)
                    if(ie .eq. nepf(ic)) then
                        iv = nloop(ie, ic)
                        i1 = nloop(1, ic)
                        !if(.not. connect(iv, i1)) then
                        !endif
                    else if(ie .eq. nepf(ic)-1) then
                        iv = nloop(ie, ic)
                        iv1 = nloop(ie+1, ic)
                        !if(.not. connect(iv, iv1)) then
                        !endif
                    else
                        iv = nloop(ie, ic)
                        iv1 = nloop(ie+1, ic)
                        if(.not. connect(iv,iv1)) then
                            do je=ie+2, nepf(ic)
                                jv = nloop(je, ic)
                                if(connect(iv, jv)) then
                                    nloop(ie+1, ic) = jv
                                    nloop(je, ic) = iv1
                                    exit
                                endif
                                !if(je .eq. nepf(ic)) then
                                !endif
                            enddo
                        endif
                    endif
                enddo
            endif
        enddo

        do ic=1, nc
            if(nepf(ic) .ne. 0) then
                do j=1, nepf(ic)
                    ivs = nloop(j, ic)
                    if(j.eq.nepf(ic))then
                        ive = nloop(1,ic)
                    else
                        ive = nloop(j+1,ic)
                    end if
                    sleng(j, ic) = sqrt((v(ivs,1)-v(ive,1))**2+(v(ivs,2)-v(ive,2))**2+(v(ivs,3)-v(ive,3))**2)
                    tleng(ic)=tleng(ic) + sleng(j,ic)
                    x1 = v(ivs,1) - v(nloop(1,ic),1)
                    y1 = v(ivs,2) - v(nloop(1,ic),2)
                    z1 = v(ivs,3) - v(nloop(1,ic),3)
                    x2 = v(ive,1) - v(nloop(1,ic),1)
                    y2 = v(ive,2) - v(nloop(1,ic),2)
                    z2 = v(ive,3) - v(nloop(1,ic),3)
                    area(ic) = area(ic) + 0.5*sqrt( (y1*z2-z1*y2)**2 + (z1*x2-z2*x1)**2 + (x1*y2-x2*y1)**2 )
                enddo
                tarea = tarea + area(ic)
                vvol(i) = vvol(i)+area(ic)*sqrt(p(ic,4))/6
            endif
        enddo
        sumvol = sumvol + vvol(i)

        !****************************************************************
        !       drop small faces / edges, optional
        !****************************************************************

        avgarea = tarea/nf
        do  ic=1, nc
            if(nepf(ic) .ne. 0) then
                if((area(ic) .ne. 0) .and.(area(ic) .lt. atol*tarea)) then
                    exit
                endif
                avglen = tleng(ic)/real(nepf(ic))
                !do j=1, nepf(ic)
                    !if((sleng(j,ic) .ne. 0.0) .and. (sleng(j,ic) .lt. tltol*avglen)) then
                    !endif
                !enddo
            endif
        enddo

        ! sccx(i) is the coordinate of the center atom
        mscce(i) = id(i)
        sccx(i) = pos(i,1)
        sccy(i) = pos(i,2)
        sccz(i) = pos(i,3)
        ! Move the center atoms according to the periodical condition
        ! nx=sccx(i)+sccx(i)
        ! sccx(i)=sccx(i)-nx
        ! ny=sccy(i)+sccy(i)
        ! sccy(i)=sccy(i)-ny
        ! nz=sccz(i)+sccz(i)
        ! sccz(i)=sccz(i)-nx

        do ic=1, nc
            if(nepf(ic) .ne. 0) then
                do j=1, nepf(ic)
                    if(sleng(j,ic).ne.0) then
                        nedges(ic,i) = nedges(ic,i)+1
                    endif
                enddo
                if(nedges(ic,i).ne.0) then
                    nnab(i) = nnab(i) + 1
                    ! by sghao, note that ic!=nnab(i) which causes a bug later on
                    ! the zero in nedges needs to be bubbled out!!!!
                    ! this bug is fixed by sghao
                    nablst(nnab(i),i) = mtag(ic)
                    ! Add by sywang revised by sghao
                    !     scx(i,j) is the coordinate of the center atom'
                    !     neighbors
                    msce(i,nnab(i))=id(mtag(ic))
                    scx(i,nnab(i))=pos(mtag(ic),1)
                    scy(i,nnab(i))=pos(mtag(ic),2)
                    scz(i,nnab(i))=pos(mtag(ic),3)

                    ! Move the atoms according to the periodical condition
                    ! For x:
                    xij=-(sccx(i)-scx(i,nnab(i)))
                    nx=xij+xij
                    scx(i,nnab(i))=scx(i,nnab(i))-nx
                    ! For y:
                    yij=-(sccy(i)-scy(i,nnab(i)))
                    ny=yij+yij
                    scy(i,nnab(i))=scy(i,nnab(i))-ny
                    ! For z:
                    zij=-(sccz(i)-scz(i,nnab(i)))
                    nz=zij+zij
                    scz(i,nnab(i))=scz(i,nnab(i))-nz
                endif
            endif
        enddo

        iicc=0
        do iic=1,nc
            if(nedges(iic,i).gt.0)then
                iicc=iicc+1
                nedges1(iicc)=nedges(iic,i)
            endif
        enddo
        do iic=1,iicc
            nedges(iic,i)=nedges1(iic)
        enddo
    enddo

    volratio = sumvol / vol
    !write(*,*) "percentages of volume counted: ", volratio
    write(*,*) "Number of bad VP = ", nbad

end subroutine vtanal


subroutine sort!(nc, p, mtag)
    use shared_data
    implicit none
    !integer, intent(in) :: nc
    !real, intent(inout) :: p(:,:)
    !integer, intent(inout) :: mtag(:)
    integer :: itag, i, j
    real :: pi1, pi2, pi3, pi4 ! Temp vars for sorting purposes.

    do i=1,nc-1
        do j=i+1,nc
            if(p(i,4).gt.p(j,4))then
                pi1 = p(i,1)
                pi2 = p(i,2)
                pi3 = p(i,3)
                pi4 = p(i,4)
                itag = mtag(i)
                p(i,1) = p(j,1)
                p(i,2) = p(j,2)
                p(i,3) = p(j,3)
                p(i,4) = p(j,4)
                mtag(i) = mtag(j)
                p(j,1) = pi1
                p(j,2) = pi2
                p(j,3) = pi3
                p(j,4) = pi4
                mtag(j) = itag
            end if
        enddo
    enddo

end subroutine sort


subroutine work!(maxver, maxepf, noi, nc, tol, p, v, nepf, nloop, mvijk, nv, nf, ne, nbad, ibad)
    use shared_data
    implicit none
    !integer, intent(in) :: nc, maxver, maxepf, noi
    !real, intent(in) :: tol
    !real, dimension(:,:), intent(in) :: p
    !integer, dimension(:), intent(inout) :: nepf, ibad
    !integer, dimension(:,:), intent(inout) :: nloop, mvijk
    !integer, intent(inout) :: nv, nf, ne, nbad
    !real, dimension(:,:), intent(inout) :: v
    integer :: i, j, k, iv, l
    real :: ai, bi, ci, di, aj, bj, cj, dj, ab, bc, ca, da, db, dc, ak, bk, ck, dk
    real :: det, detinv, vxijk, vyijk, vzijk
    logical :: ok

    if(nc.lt.4) then
        write(*,*) 'less than 4 points given to work', nc
        stop
    endif

    iv = 0
    do i=1, nc-2
        ai = p(i,1)
        bi = p(i,2)
        ci = p(i,3)
        di = -p(i,4)
        do j=i+1, nc-1
            aj =  p(j,1)
            bj =  p(j,2)
            cj =  p(j,3)
            dj = -p(j,4)
            ab = ai * bj - aj * bi
            bc = bi * cj - bj * ci
            ca = ci * aj - cj * ai
            da = di * aj - dj * ai
            db = di * bj - dj * bi
            dc = di * cj - dj * ci
            do k=j+1, nc
                ak =  p(k,1)
                bk =  p(k,2)
                ck =  p(k,3)
                dk = -p(k,4)
                det = ak * bc + bk * ca + ck * ab
                if ( abs ( det ) .gt. tol ) then
                    detinv = 1.0 / det
                    vxijk = ( - dk * bc + bk * dc - ck * db ) * detinv
                    vyijk = ( - ak * dc - dk * ca + ck * da ) * detinv
                    vzijk = (   ak * db - bk * da - dk * ab ) * detinv
                    ok = .true.
                    l=1
                    do while( ok .and. (l .le. nc) )
                        if((l.ne.i).and.(l.ne.j).and.(l.ne.k)) then
                            ok=((p(l,1)*vxijk+p(l,2)*vyijk+p(l,3)*vzijk).le.p(l,4))
                        endif
                        l=l+1
                    enddo
                    if(ok) then
                        iv = iv + 1
                        if(iv .gt. maxver) stop 'too many vertices'
                        mvijk(iv,1)  = i
                        mvijk(iv,2)  = j
                        mvijk(iv,3)  = k
                        v(iv,1) = 0.5 * vxijk
                        v(iv,2) = 0.5 * vyijk
                        v(iv,3) = 0.5 * vzijk
                    endif
                endif
            enddo
        enddo
    enddo

    nv = iv
    if(nv.lt.4) stop 'less than 4 vertices found in work'

    ! These are arrays
    nepf = 0
    nloop = 0
    !write(*,*) "Number of vertices = ", nv
    do iv=1, nv ! mvijk is set in the loop above.
        nepf(mvijk(iv,1)) = nepf(mvijk(iv,1)) + 1
        nepf(mvijk(iv,2)) = nepf(mvijk(iv,2)) + 1
        nepf(mvijk(iv,3)) = nepf(mvijk(iv,3)) + 1
        if(nepf(mvijk(iv,1)).gt.maxepf) stop 'epf>maxepf'
        if(nepf(mvijk(iv,2)).gt.maxepf) stop 'epf>maxepf'
        if(nepf(mvijk(iv,3)).gt.maxepf) stop 'epf>maxepf'
        nloop(nepf(mvijk(iv,1)),mvijk(iv,1)) = iv
        nloop(nepf(mvijk(iv,2)),mvijk(iv,2)) = iv
        nloop(nepf(mvijk(iv,3)),mvijk(iv,3)) = iv
    enddo

    nf = 0
    ne = 0
    do i=1, nc
        if(nepf(i) .gt. 0) nf = nf + 1
        ne = ne + nepf(i)
    enddo
    ! I am pretty sure a polyhedra can have an odd number of edges??? - Jason
    ! See the Futher Comments section of:
    ! http://mathcentral.uregina.ca/mp/previous2008/mar09sol.php
    ! Still, maybe we don't want this to be possible.
    ! Ahhhhh but maybe this code generates 2e, see line below if statement.
    if(mod(ne,2).ne.0)then
        write(*,*)ne
        do iv=1,nv
            write(*,*) mvijk(iv,:)
        enddo
    endif
    ne = ne/2
    if((nv-ne+nf).ne.2)then
        nbad = nbad + 1
        ibad(nbad) = noi
    endif
end subroutine work


subroutine outvt!(natoms, nsp, id, nedges, nablst, nnabsp, nnab, indx3, indx4, indx5, indx6, indx7, indx8, indx9, indx10, vvol)
    use shared_data
    implicit none
    integer, parameter :: nind=100,ncn=20,idct=2
    real, parameter :: pi=3.14159265
    !integer, intent(in) :: natoms, nsp
    !integer, dimension(:), intent(in) :: id
    !integer, dimension(:), intent(inout) :: indx3, indx4, indx5, indx6, indx7, indx8, indx9, indx10, nnab
    !integer, dimension(:,:), intent(inout) :: nnabsp
    !integer, dimension(:,:), intent(in) :: nedges, nablst
    !real, dimension(:), intent(in) :: vvol
    integer :: i, j, iind, istat!, ncluster1, ncluster2, ncluster3
    logical :: add, ins

    integer :: indlst(nind,5),icnlst(ncn)
    character(len=256) :: outbase, outfile

    outfile = "vor.out"
    outbase = ""

    outfile = trim(outbase) // '_index.out'
    !open(unit=3,file=outfile,status='unknown', access='sequential',form='formatted', iostat=istat)
    open(unit=3,file=outfile,status='replace', access='sequential',form='formatted', iostat=istat)
    if(istat .ne. 0) then
        write (*,*) "Error opening index output file, "," ",outfile, " status=",istat
        stop
    endif

    777 format(26I3)
    ! These are arrays
    indlst=0
    icnlst=0

    !ncluster1=0
    !ncluster2=0
    !ntcluster=0
    do i=1,natoms
        do j=1,nsp
            nnabsp(i,j)=0
        enddo
        indx3(i)=0
        indx4(i)=0
        indx5(i)=0
        indx6(i)=0
        indx7(i)=0
        indx8(i)=0
        indx9(i)=0
        indx10(i)=0

        do j=1,nnab(i)
            nnabsp(i,id(nablst(j,i)))=nnabsp(i,id(nablst(j,i)))+1
            if (nedges(j,i).eq.3) then
                indx3(i)=indx3(i)+1
            elseif (nedges(j,i).eq.4) then
                indx4(i)=indx4(i)+1
            elseif (nedges(j,i).eq.5) then
                indx5(i)=indx5(i)+1
            elseif (nedges(j,i).eq.6) then
                indx6(i)=indx6(i)+1
            elseif (nedges(j,i).eq.7) then
                indx7(i)=indx7(i)+1
            elseif (nedges(j,i).eq.8) then
                indx8(i)=indx8(i)+1
            elseif (nedges(j,i).eq.9) then
                indx9(i)=indx9(i)+1
            elseif (nedges(j,i).eq.10) then
                indx10(i)=indx10(i)+1
            endif
        enddo
        j=nnab(i)

        !ntcluster=ntcluster+1

        81     FORMAT(A2,I4,3F16.10)
        84     format(A2,3f16.10)
        82     FORMAT(I2,3F16.10)
        83     format(25I4)
        write(3,"(6I6,8I4,F8.4)")i,id(i),nnab(i),nnabsp(i,1),nnabsp(i,2),nnabsp(i,3),indx3(i),indx4(i),indx5(i),indx6(i),indx7(i),indx8(i),indx9(i),indx10(i),vvol(i)
        if(id(i).eq.idct)then
            icnlst(nnab(i))=icnlst(nnab(i))+1
            do iind=1,nind
                add=indlst(iind,1).eq.0.and.indlst(iind,2).eq.0.and.indlst(iind,3).eq.0.and.indlst(iind,4).eq.0
                ins=indlst(iind,1).eq.indx3(i).and.indlst(iind,2).eq.indx4(i).and.indlst(iind,3).eq.indx5(i).and.indlst(iind,4).eq.indx6(i)
                if(add)then
                    indlst(iind,1)=indx3(i)
                    indlst(iind,2)=indx4(i)
                    indlst(iind,3)=indx5(i)
                    indlst(iind,4)=indx6(i)
                    indlst(iind,5)=indlst(iind,5)+1
                    exit
                elseif(ins) then
                    indlst(iind,5)=indlst(iind,5)+1
                    exit
                endif
            enddo
        endif
    enddo
    
    close(3)

    ! There was some stuff down here, and a bit more above too, that I left out
    ! becuse it didnt do anything. Looked like debugging code.

end subroutine outvt

end module vor_hrmc
