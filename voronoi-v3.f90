! Modified to treat ternary - JWH 09/08/2009
! Converted almost completely to more modern F90-style code PMV 7/18/12
! Modified to accept command line arguments and variable file names PMV 7/18/12

program main
  implicit real(a-h,o-z)  
  parameter (n0=50000)
  parameter (nsp=3)         ! # of species
  parameter (maxcan=100)     ! max # of nearest neighbor candidate
  parameter (maxver=100)    ! max # of vertices per cell
  parameter (maxepf=20)     ! max # of edges per face
  parameter  (atol=0.01)    ! tolerance-faces to ignore, optional
  parameter  (tltol=0.01)   ! tolerance-edges to ignore, optional
  parameter  (vtol=1e-6)    ! tolerance-face intersection
  parameter  (tol=0.01)     !tolerance-faces to ignore,optional add bywsy
  character*2 atomname1,atomname2,atomname3
  common /paraatom/ ndata, n, np1, np2
  common /parasize/ a, rcut
  common /posi/ s(n0,3)     ! atomic coordinate,
  common /radii/ rsp(nsp)   ! atomic radii
  common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
  common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
  common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp) ! the species of the atoms
  common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
  common /volume/ vvol(n0)  ! vor cell volume
  common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0), indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
  common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
  common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
  common /atomtag/ mtag(maxcan) ! tag of atom after sorting
  common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
  common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
  common /aname/ atomname1,atomname2,atomname3
  common /current/ now      ! sghao to record the current NO. of snapshot      
  common /badatom/ nbad,ibad(n0)
  common /nbrsp/ nnabsp(n0,nsp) ! Number of neighbor species???
  integer,dimension(:), allocatable :: indx30,indx40,indx50,indx60, indx70,indx80,indx90,indx100,nbr10,nbr20,nbr30
  integer,dimension(:), allocatable :: indx3p,indx4p,indx5p,indx6p, indx7p,indx8p,indx9p,indx10p,nbr1p,nbr2p,nbr3p
  integer,dimension(:), allocatable :: ncount,ncount1,ncount2,ncount3,ntype

  ! filename and command line variables
  character(len=512) :: c, param_filename, model_filename, outbase, outfile
  integer :: len, istat
  
  atomname1='Zr'
  atomname2='Cu'
  atomname3='Al'
  
  ! parse command line for model name
  if(command_argument_count() /= 3) then
     write (*,*) 'vor invoked with wrong number of arguments.  Proper usage is:'
     write (*,*) ' '
     write (*,*) '     vor <parameters file>  <input model file>  <output base>'
     write (*,*) ' '
     stop
  else
     call get_command_argument(1, c, len, istat)
     if (istat == 0) then
        param_filename = trim(c)
     else
        write (*,*) 'Cannot read parameter file name.  Exiting.'
        stop
     endif
     call get_command_argument(2, c, len, istat)
     if(istat == 0) then
        model_filename = trim(c)
     else
        write (*,*) 'Cannot read model file name.  Exiting.'
        stop
     endif
     call get_command_argument(3, c, len, istat)
     if(istat == 0) then
        outbase = trim(c)
     else
        write (*,*) 'Cannot read output basename.  Exiting.'
        stop
     endif
  endif


  open(unit=2,file=model_filename, iostat=istat)
  if(istat .ne. 0) then
     write(*,*) "Error opening model file, "," ",model_filename," status=",istat
     stop
  endif
  open(unit=4,file=param_filename, iostat=istat)
  if(istat .ne. 0) then
     write(*,*) "Error opening parameters ilfe, "," ",param_filename," status=",istat
     stop
  endif

  outfile = trim(outbase) // '_index.out'
  open(unit=3,file=outfile,status='unknown', access='sequential',form='formatted', iostat=istat)  
  if(istat .ne. 0) then
     write (*,*) "Error opening index output file, "," ",outfile, " status=",istat
     stop
  endif

  outfile = trim(outbase) // '_stat.out'
  open(unit=8,file=outfile,status='replace', iostat=istat)
  if(istat .ne. 0) then
     write (*,*) "Error opening stat output file, "," ",outfile, " status=",istat
     stop
  endif
  
  read(4,*)
  read(4,*) ndata, n, np1, np2
  read(4,*)
  read(4,*) a, rcut 
  allocate(indx30(ndata*n))
  allocate(indx40(ndata*n))
  allocate(indx50(ndata*n))
  allocate(indx60(ndata*n))
  allocate(indx70(ndata*n))
  allocate(indx80(ndata*n))
  allocate(indx90(ndata*n))
  allocate(indx100(ndata*n))
  allocate(nbr10(ndata*n))
  allocate(nbr20(ndata*n))
  allocate(nbr30(ndata*n))
  allocate(indx3p(ndata*n))
  allocate(indx4p(ndata*n))
  allocate(indx5p(ndata*n))
  allocate(indx6p(ndata*n))
  allocate(indx7p(ndata*n))
  allocate(indx8p(ndata*n))
  allocate(indx9p(ndata*n))
  allocate(indx10p(ndata*n))
  allocate(nbr1p(ndata*n))
  allocate(nbr2p(ndata*n))
  allocate(nbr3p(ndata*n))
  allocate(ncount(ndata*n))
  allocate(ncount1(ndata*n))
  allocate(ncount2(ndata*n))
  allocate(ncount3(ndata*n))
  allocate(ntype(ndata*n))
  
  do jj=1,ndata
     now=jj
     nbad=0
     do kk=1,n
        ibad(kk)=0
     enddo
     !         read(2,*)
     !         read(2,*)
     do ii=1,N
        
	read(2,*) s(ii,1),s(ii,2),s(ii,3) 
        s(ii,1)=0.5*(s(ii,1)+1.0) 	!watch input model format -JWH
        s(ii,2)=0.5*(s(ii,2)+1.0)
        s(ii,3)=0.5*(s(ii,3)+1.0)
        
     enddo

80   FORMAT(2X,3F16.10)	                                        
      
     !	rsp(1)=1.58
     !  rsp(2)=1.27
     !  rsp(1)=2.16
     !  rsp(2)=1.57
     
     rsp(1)=1.0             !Hao set the radii all to zero, so I will too. - JWH   
     rsp(2)=1.0
     rsp(3)=1.0

     !	rsp(1) = 1.58		!atomic radii of zr, cu, al  - JWH 09/08/2009
     !	rsp(2) = 1.27
     !	rsp(3) = 1.43
     

     do i=1,n			!fixed for ternary - JWH 09/08/2009
        if(i.le.np1)then
           id(i)=1
        else
           if(i.gt.np1.and.i.le.np1+np2)then
              id(i)=2
           else
              id(i)=3
           endif
        endif
     enddo

     nid(1)=np1		  !total number of the first kind atom
     nid(2)=np2           !total number of the second kind atom 
     nid(3)=n-np1-np2	  !total number of the third kind atom 
     
     call initsys
     call vtanal
     call outvt
     i=0
     do ii=(jj-1)*n+1,jj*n
        i=i+1
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
        
        if(i.le.np1)then 		!fixed for ternary - JWH 09/08/2009
               ntype(i)=1
        else
               if(i.gt.np1.and.i.le.np1+np2)then
                  ntype(i)=2
               else
                  ntype(i)=3
               endif
        endif
        
     enddo
  enddo
   
     ! do the statistics
     iacu=0
     do i=1,ndata*n
        iflg=0
        do j=1,iacu
           if(indx30(i).eq.indx3p(j).and.indx40(i).eq.indx4p(j).and.indx50(i).eq.indx5p(j).and.indx60(i).eq.indx6p(j).and.indx70(i).eq.indx7p(j).and. indx80(i).eq.indx8p(j).and.indx90(i).eq.indx9p(j).and.indx100(i).eq.indx10p(j))then
              iflg=1
              ncount(j)=ncount(j)+1
              nbr1p(j)=nbr1p(j)+nbr10(i)
              nbr2p(j)=nbr2p(j)+nbr20(i)
              nbr3p(j)=nbr3p(j)+nbr30(i)
              
              if(ntype(i).eq.1)then		!fixed for ternary - JWH 09/08/2009
                 ncount1(j)=ncount1(j)+1
              else
                 if(ntype(i).eq.2)then
                    ncount2(j)=ncount2(j)+1
                 else
                    ncount3(j)=ncount3(j)+1			
                 endif
              endif
              
              !               goto 100
           endif
        enddo
        if(iflg.eq.0)then
           iacu=iacu+1
           ncount(iacu)=1
           
           if(ntype(i).eq.1)then		!fixed for ternary - JWH 09/08/2009
              ncount1(iacu)=1
           else
              if(ntype(i).eq.2)then		
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

     !fixed for ternary - JWH 09/08/2009		
     write(8,*) "No., total, center1, center2, center3, n3, n4, n5, n6, n7, n8, n9, n10, nbr-1, nbr-2, nbr-3"
     do i=1,iacu
        write(8,111) i, ncount(i), ncount1(i), ncount2(i),ncount3(i), indx3p(i),indx4p(i),indx5p(i),indx6p(i),indx7p(i),indx8p(i),indx9p(i),indx10p(i), nbr1p(i),nbr2p(i),nbr3p(i)
     enddo
111  format(5I8,8I3,3I10)
     
     
     close(2)

end program main


  !     ==================================================  
  !     initialize system and global variables
  !     ==================================================
  subroutine initsys
    implicit real(a-h,o-z) 
    parameter (n0=50000)
    parameter (nsp=3)   ! # of species
    parameter (maxcan=100)  ! max # of nearest neighbor candidate
    parameter (maxver=100) ! max # of vertices per cell
    parameter (maxepf=20) ! max # of edges per face
    parameter  (atol=0.01) ! tolerance-faces to ignore, optional
    parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
    parameter  (vtol=1e-6) ! tolerance-face intersection
    parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
    character*2 atomname1,atomname2,atomname3
    common /paraatom/ ndata, n, np1, np2
    common /parasize/ a, rcut
    common /posi/ s(n0,3) ! atomic coordinate,
    common /radii/ rsp(nsp) ! atomic radii
    common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
    common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
    common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
    common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
    common /volume/ vvol(n0) ! vor cell volume
    common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
    common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
    common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
    common /atomtag/ mtag(maxcan) ! tag of atom after sorting
    common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
    common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
    common /aname/ atomname1,atomname2,atomname3
    common /current/ now ! sghao to record the current NO. of snapshot      
    common /badatom/ nbad,ibad(n0)
    
    do i=1,n
       nnab(i)=0
       vvol(i)=0.0
       do j=1,maxcan
          nablst(j,i)=0
          nedges(j,i)=0
       end do
    end do
    
  end subroutine initsys
  
     
!     ==================================================                                
!     output results
!     ==================================================
  subroutine outvt
    implicit real(a-h,o-z)  
    parameter (n0=50000)
    parameter (nsp=3)   ! # of species
    parameter (maxcan=100)  ! max # of nearest neighbor candidate
    parameter (maxver=100) ! max # of vertices per cell
    parameter (maxepf=20) ! max # of edges per face
    parameter  (atol=0.01) ! tolerance-faces to ignore, optional
    parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
    parameter  (vtol=1e-6) ! tolerance-face intersection
    parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
    character*2 atomname1,atomname2,atomname3
    common /paraatom/ ndata, n, np1, np2
    common /parasize/ a, rcut
    common /posi/ s(n0,3) ! atomic coordinate,
    common /radii/ rsp(nsp) ! atomic radii
    common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
    common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
    common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
    common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
    common /volume/ vvol(n0) ! vor cell volume
    common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
    common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
    common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
    common /atomtag/ mtag(maxcan) ! tag of atom after sorting
    common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
    common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
    common /aname/ atomname1,atomname2,atomname3
    common /current/ now ! sghao to record the current NO. of snapshot      
    common /badatom/ nbad,ibad(n0)

    parameter(nind=100,ncn=20,pi=3.14159265,idct=2)
    parameter(list=1)
    integer ilst(5,list)
    
    integer indlst(nind,5),icnlst(ncn)
    common /nbrsp/ nnabsp(n0,nsp)
    logical add,ins,yes
        

777 format(26I3)

    do i=1,nind
       indlst(i,:)=0
    end do
    do i=1,ncn
       icnlst(i)=0
    end do
    
    ncluster1=0
    ncluster2=0
    ntcluster=0
    do i=1,n
       do j=1,nsp
          nnabsp(i,j)=0
       end do
       indx3(i)=0
       indx4(i)=0
       indx5(i)=0
       indx6(i)=0
       ! by sghao
       indx7(i)=0
       indx8(i)=0
       indx9(i)=0
       indx10(i)=0
       !!! by sghao

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

       enddo ! by sghao
       j=nnab(i)

       !       do iii=1,list
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!   Choose the cluster typy according to the Voronoi index

!      if((indx3(i).eq.0.and.indx4(i).eq.2.and.indx5(i).eq.8.and.indx6(i)
!     X.eq.2).or.(indx3(i).eq.0.and.indx4(i).eq.3.and.indx5(i).eq.6.and
!     X.indx6(i).eq.3).or.(indx3(i).eq.0.and.indx4(i).eq.0.and.indx5(i).
!     Xeq.12.and.indx6(i).eq.0)) then 
!       if(indx3(i).eq.ilst(1,iii).and.indx4(i).eq.ilst(2,iii).and.
!     x   indx5(i).eq.ilst(3,iii)
!     x  .and.indx6(i).eq.ilst(4,iii).and.nnab(i).eq.ilst(5,iii)) then
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
! here actually j=nnab(i), which is the number of ith atom

       ntcluster=ntcluster+1
      

      

81     FORMAT(A2,I4,3F16.10)	                                        	                                        
84     format(A2,3f16.10)
82     FORMAT(I2,3F16.10)     	                                        
83     format(25I4)
       write(3,"(6I6,8I3,F8.4)")i,id(i),nnab(i),nnabsp(i,1),nnabsp(i,2),nnabsp(i,3),indx3(i),indx4(i),indx5(i),indx6(i),indx7(i),indx8(i),indx9(i),indx10(i),vvol(i)			 ! added by sghao

       !       write(*,*)vvol(i)

       if(id(i).eq.idct)then
          
          icnlst(nnab(i))=icnlst(nnab(i))+1
          do iind=1,nind
             add=indlst(iind,1).eq.0.and.indlst(iind,2).eq.0 .and.indlst(iind,3).eq.0.and.indlst(iind,4).eq.0
             ins=indlst(iind,1).eq.indx3(i).and.indlst(iind,2).eq.indx4(i).and.indlst(iind,3).eq.indx5(i).and.indlst(iind,4).eq.indx6(i)
             if(add)then
                indlst(iind,1)=indx3(i)
                indlst(iind,2)=indx4(i)
                indlst(iind,3)=indx5(i)
                indlst(iind,4)=indx6(i)
                indlst(iind,5)=indlst(iind,5)+1
                exit
             elseif(ins)then
                indlst(iind,5)=indlst(iind,5)+1
                exit
             endif
          end do
       end if
    end do

    acn=0
    do i=1,ncn
       acn=acn+icnlst(i)*i
    end do
    acn=real(acn)/nid(idct)
    do i=1,ncn
       !       write(*,*)i,real(icnlst(i))/nid(idct)
    end do
    do i=1,nind
       if(indlst(i,5).ne.0)then
          !        write(*,*)indlst(i,:),real(indlst(i,5))/nid(idct)         
       endif
    end do

      

    return
  end subroutine outvt

!     ==================================================
!     voronoi analysis
!     ==================================================
  subroutine vtanal 
    implicit real(a-h,o-z)  
    common /nowatom/ noi
    logical connect
    parameter (n0=50000)
    parameter (nsp=3)   ! # of species
    parameter (maxcan=100)  ! max # of nearest neighbor candidate
    parameter (maxver=100) ! max # of vertices per cell
    parameter (maxepf=20) ! max # of edges per face
    parameter  (atol=0.01) ! tolerance-faces to ignore, optional
    parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
    parameter  (vtol=1e-6) ! tolerance-face intersection
    parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
    character*2 atomname1,atomname2,atomname3
    common /paraatom/ ndata, n, np1, np2
    common /parasize/ a, rcut
    common /posi/ s(n0,3) ! atomic coordinate,
    common /radii/ rsp(nsp) ! atomic radii
    common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
    common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
    common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
    common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
    common /volume/ vvol(n0) ! vor cell volume
    common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
    common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
    common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
    common /atomtag/ mtag(maxcan) ! tag of atom after sorting
    common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
    common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
    common /aname/ atomname1,atomname2,atomname3
    common /current/ now ! sghao to record the current NO. of snapshot      
    common /badatom/ nbad,ibad(n0)
    ! by sghao
    dimension nedges1(maxcan)
    
    sumvol=0.0
    volratio=0.0
    vol=a**3
    
    do i = 1, n
       ic=0
       noi=i
       do j = 1, n
          if(j.ne.i)then
             rxij = s(j,1)-s(i,1)
             ryij = s(j,2)-s(i,2)
             rzij = s(j,3)-s(i,3)
             
             rxij = rxij - anint ( rxij )	!The only anint used. - JWH
             ryij = ryij - anint ( ryij )
             rzij = rzij - anint ( rzij )
             
             ratio=2*rsp(id(i))/(rsp(id(i))+rsp(id(j)))       !With all radii = 0.0, the ratio is always zero. - JWH
             !        the above line implements weighted voronoi analysis
             
             rxij=rxij*ratio*a			!This is where the a is used. -JWH
             ryij=ryij*ratio*a
             rzij=rzij*ratio*a
             rijsq = rxij**2 + ryij**2 + rzij**2
             rcutsq=rcut**2
             !         rcutsq=(rcut*2*rsp(id(i))/(rsp(1)+rsp(2)))**2		 !Was ignored by Hao. - JWH 09/09/2009
             !        the above line sets cutoff for each species, optional
             if ( rijsq .lt. rcutsq ) then 
                
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
       !      candidates have been selected
       nc = ic
       write(*,*) "Number of candidates = ", nc
       !      sort into ascending order of distance
       call sort
       call work
       !      perform voronoi analysis
!****************************************************************
!             sort vertices ring in faces and 
!          calculate edge length and face area
!****************************************************************
       tarea = 0.0
       avglen = 0.0
       avgarea = 0.0
       do ic = 1, maxcan
          do jedge = 1, maxepf
             sleng(jedge,ic) = 0.0
          end do
          tleng(ic) = 0.0
          area(ic) = 0.0
       end do
          
       do ic = 1, nc
          !if(nepf(ic) .eq. 1) write(*,*) "Careful! There could be an error here!!!" ! Jason Maldonis 11/14/2013
          !if (nepf(ic) .ne. 0 .and. nepf(ic) .ne. 1) then  ! JASON MALDONIS 11/14/2013 Changed from 'if (nepf(ic) .ne. 0) then'
          if (nepf(ic) .ne. 0) then
             do ie = 1, nepf(ic)
                if(ie.eq.nepf(ic))then
                   iv=nloop(ie,ic)
                   i1=nloop(1,ic)
    if(iv .eq. i1) write(*,*) "WARNING ERROR?", ie, ic, nepf(ic)
    if(iv .eq. i1) write(*,*) nepf, nc
                   if(.not.connect(iv,i1))then
                      !           write(27,*) 'not a loop-0'
                      !           write(27,*) now, noi
                   endif
                else if(ie.eq.nepf(ic)-1)then
                   iv=nloop(ie,ic)
                   iv1=nloop(ie+1,ic)
    if(iv .eq. iv1) write(*,*) "WARNING ERROR?"
                   if(.not.connect(iv,iv1))then
                      !           write(27,*) 'not a loop-1'
                      !           write(27,*) now,noi
                   endif
                else
                   iv=nloop(ie,ic)
                   iv1=nloop(ie+1,ic)
    if(iv .eq. iv1) write(*,*) "WARNING ERROR?"
                   if(.not.connect(iv,iv1))then
                      do je=ie+2,nepf(ic)
                         jv=nloop(je,ic)
    if(iv .eq. jv) write(*,*) "WARNING ERROR?"
                         if(connect(iv,jv))then
                            nloop(ie+1,ic)=jv
                            nloop(je,ic)=iv1
                            exit
                         end if
                         if(je.eq.nepf(ic))then
                            !              write(27,*) 'not a loop-2'
                            !              write(27,*) now, noi
                         endif
                      end do
                   end if
                end if
             end do
          end if
       enddo
             
       do ic = 1, nc
          if (nepf(ic) .ne. 0) then  
             do j = 1, nepf(ic)
                ivs = nloop(j,ic)
                if(j.eq.nepf(ic))then
                   ive = nloop(1,ic)
                else
                   ive = nloop(j+1,ic)
                end if
                sqlen = (v(ivs,1)-v(ive,1))**2+(v(ivs,2)-v(ive,2))**2+(v(ivs,3)-v(ive,3))**2
                sleng(j, ic) = sqrt(sqlen)
                tleng(ic)=tleng(ic) + sleng(j,ic)
                x1 = v(ivs,1) - v(nloop(1,ic),1)
                y1 = v(ivs,2) - v(nloop(1,ic),2)
                z1 = v(ivs,3) - v(nloop(1,ic),3)
                x2 = v(ive,1) - v(nloop(1,ic),1)
                y2 = v(ive,2) - v(nloop(1,ic),2)
                z2 = v(ive,3) - v(nloop(1,ic),3)
                area(ic) = area(ic) + .5*sqrt((y1*z2-z1*y2)**2+(z1*x2-z2*x1)**2+(x1*y2-x2*y1)**2)
             enddo
             tarea = tarea + area(ic)
             vvol(i)=vvol(i)+area(ic)*sqrt(p(ic,4))/6
          end if
       enddo
       sumvol = sumvol + vvol(i)
       
       !****************************************************************
!             drop small faces / edges, optional
       !****************************************************************
       
       avgarea = tarea/nf
       do  ic = 1, nc
          if (nepf(ic) .ne. 0) then
             if ((area(ic) .ne. 0) .and.(area(ic) .lt. atol*tarea)) then
                !           sleng(:,ic)=0
                !          pause "droped a face!"
                exit
             end if
             avglen = tleng(ic)/real(nepf(ic))
             do j = 1, nepf(ic)
                if ((sleng(j,ic) .ne. 0.0) .and. (sleng(j,ic) .lt. tltol*avglen)) then
                   !              sleng(j,ic)=0
                   !             pause "droped an edge!"
                end if
             end do
          end if
       enddo
       
       !     sccx(i) is the coordinate of the center  atom    
       mscce(i) = id(i)	
       sccx(i)=s(i,1)
       sccy(i)=s(i,2)
       sccz(i)=s(i,3)
       ! move the center atoms according to the periodical condition 
       !          nx=sccx(i)+sccx(i)
       !          sccx(i)=sccx(i)-nx
       
       !          ny=sccy(i)+sccy(i)
       !          sccy(i)=sccy(i)-ny
       
       !          nz=sccz(i)+sccz(i)
       !          sccz(i)=sccz(i)-nx
       
       
       do ic = 1, nc
          if (nepf(ic) .ne. 0) then
             do j=1,nepf(ic)
                if(sleng(j,ic).ne.0)then
                   nedges(ic,i)=nedges(ic,i)+1
                end if
             end do
             if(nedges(ic,i).ne.0)then
                nnab(i) = nnab(i) + 1
                ! by sghao, note that ic!=nnab(i) which causes a bug later on
                ! the zero in nedges needs to be bubbled out!!!!
                ! this bug is fixed by sghao 
                !         write(17,*) i,ic,nnab(i),nedges(ic,i) 
                ! by sghao
                nablst(nnab(i),i) = mtag(ic)
                ! Add by sywang revised by sghao
                !     scx(i,j) is the coordinate of the center atom' neighbors    
                msce(i,nnab(i))=id(mtag(ic))
                scx(i,nnab(i))=s(mtag(ic),1)
                scy(i,nnab(i))=s(mtag(ic),2)
                scz(i,nnab(i))=s(mtag(ic),3) 
                
                ! move the atoms according to the periodical condition
                !  for x
                xij=-(sccx(i)-scx(i,nnab(i)))
                nx=xij+xij
                scx(i,nnab(i))=scx(i,nnab(i))-nx
                !  for y   	  
                yij=-(sccy(i)-scy(i,nnab(i)))
                ny=yij+yij
                scy(i,nnab(i))=scy(i,nnab(i))-ny
                
                !  for z	        
                zij=-(sccz(i)-scz(i,nnab(i)))
                nz=zij+zij
                scz(i,nnab(i))=scz(i,nnab(i))-nz
                
                ! Add by sywang revised by sghao
                
             end if
          endif
       enddo
       
       ! added by sghao
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
       !ccccccby sghao      
       
    end do
    
!    *******************************************************************
!    ** main loop ends                                                **
!    *******************************************************************
    volratio = sumvol / vol
    !      write(*,*) "percentages of volume counted: ", volratio
                         
    return
  end subroutine vtanal
                       

!     ==================================================

!     ==================================================
   function connect(i,j)
     
     implicit real(a-h,o-z) 
     parameter (n0=50000)
     parameter (nsp=3)   ! # of species
     parameter (maxcan=100)  ! max # of nearest neighbor candidate
     parameter (maxver=100) ! max # of vertices per cell
     parameter (maxepf=20) ! max # of edges per face
     parameter  (atol=0.01) ! tolerance-faces to ignore, optional
     parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
     parameter  (vtol=1e-6) ! tolerance-face intersection
     parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
     character*2 atomname1,atomname2,atomname3
     common /paraatom/ ndata, n, np1, np2
     common /parasize/ a, rcut
     common /posi/ s(n0,3) ! atomic coordinate,
     common /radii/ rsp(nsp) ! atomic radii
     common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
     common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
     common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
     common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
     common /volume/ vvol(n0) ! vor cell volume
     common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
     common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
     common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
     common /atomtag/ mtag(maxcan) ! tag of atom after sorting
     common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
     common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
     common /aname/ atomname1,atomname2,atomname3
     common /current/ now ! sghao to record the current NO. of snapshot      
     common /badatom/ nbad,ibad(n0)
     logical connect 
     integer i,j
     np=0
     do ii=1,3
        do jj=1,3
!write(*,*) mvijk(i,ii), mvijk(j,jj)
           if(mvijk(i,ii).eq.mvijk(j,jj))then
              np=np+1
           end if
        end do
     end do
     select case (np)
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
!     ==================================================

!     construct voronoi cells by calculating vertices
!     ==================================================
   subroutine work
     !      use global 
     implicit real(a-h,o-z)  
     parameter (n0=50000)
     parameter (nsp=3)   ! # of species
     parameter (maxcan=100)  ! max # of nearest neighbor candidate
     parameter (maxver=100) ! max # of vertices per cell
     parameter (maxepf=20) ! max # of edges per face
     parameter  (atol=0.01) ! tolerance-faces to ignore, optional
     parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
     parameter  (vtol=1e-6) ! tolerance-face intersection
     parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
     character*2 atomname1,atomname2,atomname3
     common /paraatom/ ndata, n, np1, np2
     common /parasize/ a, rcut
     common /posi/ s(n0,3) ! atomic coordinate,
     common /radii/ rsp(nsp) ! atomic radii
     common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
     common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
     common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
     common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
     common /volume/ vvol(n0) ! vor cell volume
     common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
     common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
     common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
     common /atomtag/ mtag(maxcan) ! tag of atom after sorting
     common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
     common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
     common /aname/ atomname1,atomname2,atomname3
     common /current/ now ! sghao to record the current NO. of snapshot      
     common /badatom/ nbad,ibad(n0)
     common /nowatom/ noi
     logical ok
 
     if(nc.lt.4)stop'less than 4 points given to work'
     nc1 = nc - 1
     nc2 = nc - 2
     iv = 0
     do i = 1, nc2
        ai =  p(i,1)
        bi =  p(i,2)
        ci =  p(i,3)
        di = -p(i,4)
        do j=i + 1, nc1
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
           do k = j + 1, nc
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
                 l  = 1
100              if ( ok .and. ( l .le. nc ) ) then
                    if ((l.ne.i).and.(l.ne.j).and.(l.ne.k)) then
                       ok=((p(l,1)*vxijk+p(l,2)*vyijk+p(l,3)*vzijk).le.p(l,4))
                    endif
                    l = l + 1
                    goto 100
                 endif
                 if ( ok ) then
                    iv = iv + 1
                    if ( iv .gt. maxver ) stop 'too many vertices'
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
     if(nv .eq. 1) write(*,*) "WARNING!!! nv = 1!"
     ! sghao
     if(nv.lt.4)stop'less than 4 vertices found in work'
     do  i = 1, maxcan
        nepf(i) = 0
        do j=1, maxepf
           nloop(j,i)=0
        end do
     enddo
write(*,*) "Number of vertices =", nv, maxepf
     do iv = 1, nv
        nepf(mvijk(iv,1)) = nepf(mvijk(iv,1)) + 1
        if(nepf(mvijk(iv,1)).gt.maxepf)stop'epf>maxepf'
        nloop(nepf(mvijk(iv,1)),mvijk(iv,1))=iv
        nepf(mvijk(iv,2)) = nepf(mvijk(iv,2)) + 1
        if(nepf(mvijk(iv,2)).gt.maxepf)stop'epf>maxepf'
        nloop(nepf(mvijk(iv,2)),mvijk(iv,2))=iv
        nepf(mvijk(iv,3)) = nepf(mvijk(iv,3)) + 1
        if(nepf(mvijk(iv,3)).gt.maxepf)stop'epf>maxepf'
        nloop(nepf(mvijk(iv,3)),mvijk(iv,3))=iv
     enddo
do iv = 1, 75
    if(nepf(iv) .eq. 1) write(*,*) "WARNING ERROR??? in nepf", iv
enddo
     nf = 0
     ne = 0
     do i = 1, nc
        if ( nepf(i) .gt. 0 ) nf = nf + 1
        ne = ne + nepf(i)
     end do
     if(mod(ne,2).ne.0)then
        write(*,*)ne
        do iv=1,nv
           write(*,*)mvijk(iv,:)
        end do
     end if
     ne = ne / 2
     !      write(7,*) nv,ne,nf,nv-ne+nf 
     if((nv-ne+nf).ne.2)then 
        !       write(27,*)'euler error: degeneracy'
        !       write(27,*) now, noi
        nbad=nbad+1
        ibad(nbad)=noi
     endif
                       
     return
   end subroutine work
!     ==================================================

!     ==================================================
   subroutine sort 
  
     implicit real(a-h,o-z)  
     parameter (n0=50000)
     parameter (nsp=3)   ! # of species
     parameter (maxcan=100)  ! max # of nearest neighbor candidate
     parameter (maxver=100) ! max # of vertices per cell
     parameter (maxepf=20) ! max # of edges per face
     parameter  (atol=0.01) ! tolerance-faces to ignore, optional
     parameter  (tltol=0.01) ! tolerance-edges to ignore, optional
     parameter  (vtol=1e-6) ! tolerance-face intersection
     parameter  (tol=0.01)  !tolerance-faces to ignore,optional add bywsy
     character*2 atomname1,atomname2,atomname3
     common /paraatom/ ndata, n, np1, np2
     common /parasize/ a, rcut
     common /posi/ s(n0,3) ! atomic coordinate,
     common /radii/ rsp(nsp) ! atomic radii
     common /clstposi/ scx(n0,25),scy(n0,25),scz(n0,25) ! atomic coordinate of the cluster
     common /centposi/ sccx(n0),sccy(n0),sccz(n0) ! atomic coordinate of thecenter atom at cluster
     common /specatom/ mscce(n0),msce(n0,25),id(n0),nid(nsp)! the species of the atoms
     common /nghbor/ nnab(n0),nablst(maxcan,n0),nedges(maxcan,n0) ! neighbor list
     common /volume/ vvol(n0) ! vor cell volume
     common /vrnindex/ indx3(n0),indx4(n0),indx5(n0),indx6(n0),indx7(n0),indx8(n0),indx9(n0),indx10(n0) ! Voronoi indices
     common /cvef/ nc,nv,ne,nf,mvijk(maxver,3) ! # of candidates,vertices,edges,faces
     common /vertices/ p(maxcan,4),v(maxver,3) ! vertices
     common /atomtag/ mtag(maxcan) ! tag of atom after sorting
     common /concheck/ nloop(maxepf,maxcan),nepf(maxcan) ! check connection
     common /arealen/ area(maxcan),sleng(maxepf,maxcan),tleng(maxcan) ! area and length
     common /aname/ atomname1,atomname2,atomname3
     common /current/ now ! sghao to record the current NO. of snapshot      
     common /badatom/ nbad,ibad(n0)
     
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
     return
   end subroutine sort
!     ==================================================

