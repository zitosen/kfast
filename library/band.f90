! by zito 2020.12.18
      PROGRAM band_data_process
      IMPLICIT NONE
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: kpoints
      REAL(KIND=8) :: flag,flag_step,step
! 'ncol' is number of lines of each block
! 'nblock' is number of blocks
      INTEGER, PARAMETER :: ncol=299, nblock=95
      CHARACTER(100) :: temp
      INTEGER :: i,k,j
!
      ALLOCATE(kpoints(1:ncol))
      OPEN(unit=77,file='band.dat',status='old')
      DO i=1,ncol
        READ(77,*) kpoints(i)
      ENDDO
      CLOSE(77)
!
      step=0.d0
      DO WHILE(.TRUE.)
        k=0
        DO i=1,ncol
          IF(i>1) step=kpoints(i)-kpoints(i-1)
          IF((kpoints(i).GT.1.3d0) .AND. (step.GT.0.02)) THEN
            flag=kpoints(i)
            flag_step=step
            EXIT
          ENDIF
          k=i
        ENDDO
        IF(k==ncol) EXIT
        DO i=1,ncol
          IF( kpoints(i).GT.flag-0.0001 ) THEN
            kpoints(i)=kpoints(i)-flag_step
          ENDIF
        ENDDO
      ENDDO
!
      OPEN(unit=88,file='process.dat',status='unknown')
      DO j=1,nblock  
        DO i=1,ncol
          WRITE(88,"(3F15.7)") kpoints(i)
        ENDDO
        WRITE(88,*) ""
      ENDDO
      CLOSE(88)
      DEALLOCATE(kpoints)
!      
      END
