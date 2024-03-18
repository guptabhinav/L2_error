PROGRAM l2_error_pbtass

IMPLICIT NONE
!Define global variables
  INTEGER :: n1, n2
  REAL*8  :: e_cut,l2error
  REAL*8, POINTER :: fes_remd(:,:), fes_tass(:,:),error(:,:),grid(:,:)
  LOGICAL, POINTER :: ind(:,:)
 
  CALL read_inputs
  CALL align_surfaces

  CALL get_L2error(l2error)
  PRINT *, 'L2 error= ', l2error
  CALL print_error

CONTAINS

  SUBROUTINE align_surfaces
  IMPLICIT NONE

  REAL*8 :: minf_remd,minf_tass,avfe

  minf_remd=get_minfe(fes_remd)
  PRINT *, 'Minimum of REMD free energy =', minf_remd
  fes_remd(:,:)=fes_remd(:,:)-minf_remd

  minf_tass=get_minfe(fes_tass)
  PRINT *, 'Minimum of TASS free energy =', minf_tass
  fes_tass(:,:)=fes_tass(:,:)-minf_tass

  CALL get_indices(e_cut,fes_remd)

! Correct the free energy surface 
  avfe=get_average(fes_remd)
  fes_remd(:,:)=fes_remd(:,:)-avfe

  avfe=get_average(fes_tass)
  fes_tass(:,:)=fes_tass(:,:)-avfe

  END SUBROUTINE align_surfaces

!-------------------------------------------------!
  SUBROUTINE print_error
  IMPLICIT NONE
  INTEGER :: ii, i, j
!
  OPEN(101,file='l2error.dat')
!
  ii=0
  DO i=1,n1
   DO j=1,n2
     ii=ii+1
     WRITE(101,'(3F16.6)')grid(ii,1:2),error(i,j) 
   END DO
   WRITE(101,*)
  END DO
  CLOSE(101)
!
  END SUBROUTINE print_error

!-------------------------------------------------!
  SUBROUTINE get_L2error(l2error)
  IMPLICIT NONE
  REAL*8 :: l2error
!
  INTEGER :: i,j,ii
!
  ALLOCATE(error(n1,n2))
  error(:,:)=999.d0
  
  l2error=0.d0
  ii=0
  DO i=1,n1
   DO j=1,n2
    IF(ind(i,j))THEN
      ii=ii+1
      l2error=l2error+(fes_tass(i,j)-fes_remd(i,j))**2
      error(i,j)=DSQRT((fes_tass(i,j)-fes_remd(i,j))**2)
    END IF
   END DO
  END DO
  l2error=DSQRT((l2error)/DBLE(ii))
!
  END SUBROUTINE get_L2error


!-------------------------------------------------!
  SUBROUTINE get_indices(e_cut,fes)
  IMPLICIT NONE
!
  REAL*8 :: e_cut
  REAL*8 :: fes(n1,n2)
!
  INTEGER :: i, j, ii
!
  ALLOCATE(ind(n1,n2))
  ind(:,:)=.FALSE.
  ii=0
  DO i=3,n1-2
   DO j=3,n2-2 !Ignoring the boundaries
    IF(fes(i,j)<e_cut)THEN
       ii=ii+1
       ind(i,j)=.TRUE.
    END IF
   END DO
  END DO
  PRINT *, 'Number of points considered = ', ii
  END SUBROUTINE get_indices
!-------------------------------------------------!
  FUNCTION get_average(fes) RESULT(av)
  IMPLICIT NONE
  REAL*8 :: fes(n1,n2)
  INTEGER :: i,j,ii
  REAL*8 :: av
!
  av=0.d0
  ii=0
  DO i=1,n1
    DO j=1,n2
      if(ind(i,j))THEN
       ii=ii+1
       av=av+fes(i,j)
      end if
    END DO
  END DO
  av=av/DBLE(ii)
  print *, 'average fe =', av, ' (# points =',ii,' )'
!
  END FUNCTION get_average
!-------------------------------------------------!
  FUNCTION get_minfe(fes) RESULT(minf)
  IMPLICIT NONE
  REAL*8 :: fes(n1,n2)
  REAL*8 :: minf
  INTEGER :: i, j

  minf=fes(1,1)
  DO i=1,n1
    DO j=1,n2
      IF(minf>fes(i,j))minf=fes(i,j)
    END DO
  END DO
  END FUNCTION get_minfe
!-------------------------------------------------!

  SUBROUTINE read_inputs
  
  IMPLICIT NONE
  character(len=100):: file1,file2
  INTEGER :: i, j, ii
  

  WRITE(*,*)'Enter the name of the REMD free energy file'
!  READ(*,*)file1
   file1="Free_energy_remd"
  WRITE(*,*)'Enter the name of the TASS free energy file'
!  READ(*,*)file2
   file2="Free_energy_wrt_0"
  OPEN(1,file=file1)
  OPEN(2,file=file2)


  WRITE(*,*)'Enter the grid dimension along cv1 and cv2' 
  WRITE(*,*)'NOTE: Dimension should be the same for ',trim(file1), ' & ', trim(file2)

!  READ(*,*)n1,n2
  n1=100
  n2=100
  ALLOCATE(fes_remd(n1,n2))
  ALLOCATE(fes_tass(n1,n2))
  ALLOCATE(grid(n1*n2,2))
  
 ii=0
  DO i=1,n1
   DO j=1,n2
     ii=ii+1
     READ(1,*)grid(ii,1:2),fes_remd(i,j)
     READ(2,*)grid(ii,1:2),fes_tass(i,j)
   END DO
   READ(1,*) !Skipping the space
   READ(2,*) !Skipping the space
  END DO
  PRINT *, ' Finished reading ', trim(file1), ' & ', trim(file2)

  PRINT *, ' Enter the free energy cutoff'
  READ(*,*)e_cut

  END SUBROUTINE read_inputs

END PROGRAM l2_error_pbtass
