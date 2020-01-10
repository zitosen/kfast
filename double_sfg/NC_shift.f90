!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     NC_shift: calculate the shifts of Normal Coordinates              !
!               For the displaced oscillor approximation                ! 
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE NC_shift
      IMPLICIT NONE
!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Delta_Q(natom,dq)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(OUT) :: dq
      INTEGER :: i,j,k,itemp,lblock
      INTEGER :: mass(natom)
      REAL(KIND=8),DIMENSION(3,natom) :: xcoor,xcoorp
      REAL(KIND=8),DIMENSION(natom*3) :: dx
      REAL(KIND=8),DIMENSION(natom*3-6) :: redmass
      REAL(KIND=8),DIMENSION(natom*3,natom*3-6) :: lmat
      CHARACTER(100) :: filename,buffer
      LOGICAL :: ok
!
      filename="NC_shift.in"
      redmass(:)=0.d0
      dq(:)=0.d0
      INQUIRE(file=filename,exist=ok)
      IF(ok) THEN
        OPEN(unit=79,file=filename,form='formatted',status='old')
! read reduced mass
        READ(79,*) buffer
        DO i=1,natom*3-6
          READ(79,*) redmass(i)   ! in AMU 
        END DO
! read geometry of S0 state in Cartesian coordinates
        READ(79,*) buffer
        DO i=1,natom
          READ(79,*) itemp,itemp,xcoor(1:3,i)
        END DO
! read geometry of S1 state in Cartesian coordinates
        READ(79,*) buffer
        DO i=1,natom
          READ(79,*) itemp,itemp,xcoorp(1:3,i)
        END DO
! read l matrix of S0 state in Gaussian form
        READ(79,*) buffer
        IF(MOD(natom*3-6,5).NE.0) THEN
          lblock=(natom*3-6)/5+1
          DO i=1,lblock-1
            DO k=1,8 
              READ(79,*) buffer
            END DO
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,(i*5-4):(i*5))
            END DO
          END DO
          DO k=1,8   ! read the last block of l matrix
            READ(79,*) buffer
          END DO
          IF(MOD(natom*3-6,5)==1) THEN
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,natom*3-6)
            ENDDO
          ELSE IF(MOD(natom*3-6,5)==2) THEN
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,natom*3-7), &
                  lmat(j,natom*3-6)
            ENDDO
          ELSE IF(MOD(natom*3-6,5)==3) THEN
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,natom*3-8), &
                  lmat(j,natom*3-7),lmat(j,natom*3-6)
            ENDDO
          ELSE IF(MOD(natom*3-6,5)==4) THEN
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,natom*3-9), &
                  lmat(j,natom*3-8),lmat(j,natom*3-7),lmat(j,natom*3-6)
            ENDDO
          ENDIF
        ELSE
          lblock=(natom*3-6)/5
          DO i=1,lblock
            DO k=1,8
              READ(79,*) buffer
            END DO
            DO j=1,natom*3
              READ(79,*) itemp,itemp,itemp,lmat(j,(i*5-4):(i*5))
            END DO
          END DO
        END IF
!      
      ELSE
        WRITE(*,*) "File ",'"',TRIM(filename),'"', " doesn't exist."
      END IF
      CLOSE(79)
!
! calculate displacement between S0 and S1 states in Cartesian coordinates
!     
      dx=0.d0
      DO i=1,natom
        dx(i*3-2)=xcoorp(1,i)-xcoor(1,i)
        dx(i*3-1)=xcoorp(2,i)-xcoor(2,i)
        dx(i*3)  =xcoorp(3,i)-xcoor(3,i)
      ENDDO
!
! calculate normal coordinate displacements dq
!
      dq=0.d0
      DO i=1,natom*3-6
        DO j=1,natom*3
          dq(i)=dq(i)+lmat(j,i)*dx(j)
        END DO
          dq(i)=dq(i)*DSQRT(redmass(i))    ! unit: AMU^(1/2)*Angstrom
      END DO
!
      RETURN
      END SUBROUTINE Delta_Q
      END MODULE NC_shift

