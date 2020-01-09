!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     NC_shift: calculate the shifts of Normal Coordinates              !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE NC_shift
      USE contants
      USE atom_mass
      IMPLICIT NONE
!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Delta_Q(natom,freq,dq)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(OUT) :: freq,dq
      INTEGER :: i,j,itemp
      INTEGER :: mass(natom)
      REAL(KIND=8),DIMENSION(3,natom*3-6) :: xcoor,xcoorp,dx
      REAL(KIND=8),DIMENSION(natom*3-6) :: dqp,redmass,f_norm
      REAL(KIND=8),DIMENSION(natom*3-6,natom*3) :: lmatrix,lmw
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
          READ(79,*) redmass(i) 
        END DO
! read geometry of S0 state in Cartesian coordinates
        READ(79,*) buffer
        DO i=1,natom
          READ(79,*) itemp,itemp,xcoor(1,i),xcoor(2,i),xcoor(3,i)
        END DO
! read geometry of S1 state in Cartesian coordinates
        READ(79,*) buffer
        DO i=1,natom
          READ(79,*) itemp,itemp,xcoorp(1,i),xcoorp(2,i),xcoorp(3,i)
        END DO
      ELSE
        WRITE(*,*) "File ",'"',TRIM(filename),'"', " doesn't exist."
      END IF
      CLOSE(79)
!
! calculate displacement between S0 and S1 states
      dx(:,:)=xcoorp(:,:)-xcoor(:,:)
!
! the normal modes calculated by Gaussian, i.e., l matrix are read in
      filename="acetone.norm"
      INQUIRE(file=filename,exist=ok)
      IF(ok) THEN
        OPEN(unit=79,file=filename,form='formatted',status='old')
        i=1;io=0
        DO WHILE((i*5)<=(natom*3-6))
          DO WHILE(.TRUE.)
            READ(unit=79,fmt='(A100)',iostat=io) buffer
            buffer=ADJUSTL(buffer)
            IF(buffer(1:11)==nu_flag) THEN
              BACKSPACE(79)
              READ(79,*) segment,segment,freq((i*5-4):i*5)         ! read frequency
              READ(79,*) segment,segment,segment,redmass((i*5-4):i*5) ! read reduced mass
              DO j=1,3
                READ(79,'(A100)') buffer
              END DO
              DO j=1,natom*3
                READ(unit=79,fmt=*,iostat=io) segment,segment,segment, &
     &          lmatrix((i*5-4):i*5,j)
              END DO
            IF(.TRUE.) EXIT
            END IF
          END DO
          i=i+1
        END DO
        IF((i-1)*5<(natom*3-6)) THEN
          READ(unit=79,fmt='(A100)',iostat=io) buffer
          READ(unit=79,fmt='(A100)',iostat=io) buffer
          READ(79,*) segment,segment,freq(((i-1)*5+1):(natom*3-6)) ! read frequency
          READ(79,*) segment,segment,segment, &
     &               redmass(((i-1)*5+1):(natom*3-6))              ! read reduced mass
          DO j=1,3
            READ(79,'(A100)') buffer
          END DO
          DO j=1,natom*3
            READ(unit=79,fmt=*,iostat=io) segment,segment,segment, &
     &      lmatrix(((i-1)*5+1):(natom*3-6),j)
          END DO
        END IF
      ELSE
        WRITE(*,*) "File ",'"',TRIM(filename),'"', " doesn't exist."
      END IF
      CLOSE(79)
!
! calculate normal coordinate displacements: dqp !!!!!!!!!!!!!!!!!!!!!!
!
      dq=0.d0  ! unit: AMU^(1/2)*Angstrom
      dqp=0.d0 ! unit: Angstrom
      DO i=1,natom*3-6
        DO j=1,natom*3
          dqp(i)=dqp(i)+lmatrix(i,j)*dx(1,j)
        END DO
          dq(i)=dqp(i)*DSQRT(redmass(i))    ! unit: AMU^(1/2)*Angstrom
      END DO
!
      RETURN
      END SUBROUTINE Delta_Q
      END MODULE NC_shift

