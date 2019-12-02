! Summation of site-projected DOS
! 2018.9.5, by zito

      PROGRAM PDOS_SUM
      IMPLICIT NONE
      INTEGER, PARAMETER       :: sr=4,sc=8,dr=8,dc=16
      INTEGER                  :: ngrid,istart,iend,i,j
      CHARACTER(LEN=4)         :: label
      CHARACTER(LEN=20)       :: filename,buffer
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:) &
                               :: energy,s_up,s_down,p_up,p_down,d_up, &
                                  d_down
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) &
                               :: sum_dos
      REAL(KIND=dr),DIMENSION(7) &
                               :: temp
!
      WRITE(*,*) "Number of grids: "
      READ(*,*) ngrid
      WRITE(*,*) "The first and the last atom label: istart, iend"
      READ(*,*) istart, iend
      WRITE(*,*) "Summation is from ", istart, " - ", iend
      ALLOCATE(energy(1:ngrid))
      ALLOCATE(s_up(1:ngrid))
      ALLOCATE(s_down(1:ngrid))
      ALLOCATE(p_up(1:ngrid))
      ALLOCATE(p_down(1:ngrid))
      ALLOCATE(d_up(1:ngrid))
      ALLOCATE(d_down(1:ngrid))
      ALLOCATE(sum_dos(1:6,1:ngrid))
      energy=0.d0
      s_up=0.d0
      s_down=0.d0
      p_up=0.d0
      p_down=0.d0
      d_up=0.d0
      d_down=0.d0
      sum_dos=0.d0
!
      DO i=istart, iend
        WRITE(label,"(I3)") i
        filename="DOS"//TRIM(ADJUSTL(label))
        WRITE(*,*) filename, "Reading..."
!
        OPEN(UNIT=77,FILE=filename,STATUS="OLD")
        IF(i.LT.9) READ(77,*) buffer ! note the DOS1...DOS8 files
        DO j=1,ngrid
        READ(77,*) energy(j),s_up(j),s_down(j),p_up(j),p_down(j), &
                   d_up(j),d_down(j)
        ENDDO
        CLOSE(77)
!
        sum_dos(1,:)=sum_dos(1,:)+s_up(:)
        sum_dos(2,:)=sum_dos(2,:)+s_down(:)
        sum_dos(3,:)=sum_dos(3,:)+p_up(:)
        sum_dos(4,:)=sum_dos(4,:)+p_down(:)
        sum_dos(5,:)=sum_dos(5,:)+d_up(:)
        sum_dos(6,:)=sum_dos(6,:)+d_down(:)
      ENDDO
      OPEN(UNIT=7,FILE="sum_dos.dat",STATUS="UNKNOWN")
      DO i=1,ngrid
        WRITE(7,"(7F14.8)") energy(i),sum_dos(1,i),sum_dos(2,i), &
                            sum_dos(3,i),sum_dos(4,i),sum_dos(5,i),&
                            sum_dos(6,i)
      ENDDO
      CLOSE(7)
!
      DEALLOCATE(energy)
      DEALLOCATE(s_up)
      DEALLOCATE(s_down)
      DEALLOCATE(p_up)
      DEALLOCATE(p_down)
      DEALLOCATE(d_up)
      DEALLOCATE(d_down)
      END
