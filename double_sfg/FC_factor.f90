!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     atom_mass: assignment of atomic mass                              !
!                                                                       !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE atom_mass
      IMPLICIT NONE
!
      CONTAINS
      FUNCTION n2m(iso,num)
      IMPLICIT NONE
      REAL(KIND=8) :: amu(0:2,1:10),n2m
      INTEGER,INTENT(IN) :: iso,num !!iso: isotope type!
!
      n2m=0.d0
      amu(0,1:10)=(/1.007825d0,4.00260d0,7.016003d0,9.012182d0,&
     &11.009305d0,12.0d0,14.003074d0,15.994915d0,18.9984032d0,&
     &19.992435d0/)
      IF(0<iso<=2.AND.0<num<=10) THEN
        n2m=amu(iso,num)
      ELSE
        WRITE(*,'(A30)') 'Check the "atom_mass" module!'
      END IF
      END FUNCTION n2m
      END MODULE atom_mass
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     Delta_NMC: Delta Normal Mode Coordinate                           !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE Delta_NMC
      USE global
      USE atom_mass
      IMPLICIT NONE
!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Delta_Q(natom,freq,dq)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(OUT) :: freq,dq
      INTEGER :: i,j,itemp,istate,nstate,io
      INTEGER :: mass(natom)
      REAL(KIND=8) :: xcoor(3,natom*3)
      REAL(KIND=8),DIMENSION(natom*3-6) :: dqp,redmass,f_norm
      REAL(KIND=8),DIMENSION(natom*3-6,natom*3) :: lmatrix,lmw
      REAL(KIND=8),ALLOCATABLE :: dx(:,:)
      CHARACTER(100) :: filename,buffer
      CHARACTER(40) :: segment
      CHARACTER(20), PARAMETER :: begin_flag="BEGIN_", &
     &nu_flag="Frequencies"
      LOGICAL :: ok
!
      nstate=0
      filename="acetone.coor"
      INQUIRE(file=filename,exist=ok)
      IF(ok) THEN
        OPEN(unit=79,file=filename,form='formatted',status='old')
        DO WHILE(.TRUE.)
          READ(unit=79,fmt='(A100)',iostat=io) buffer
          buffer=ADJUSTL(buffer)
          IF(buffer(1:6)==begin_flag) THEN
            itemp=INDEX(buffer,"_")
            READ(buffer(itemp+1:itemp+2),'(I1)') istate
            IF(istate>nstate) nstate=istate
            DO i=1,natom
              READ(79,*) mass(i),xcoor(istate,i*3-2),xcoor(istate,i*3-1)&
     &        ,xcoor(istate,i*3)
            END DO
          END IF
          IF(io/=0) EXIT
        END DO
      ELSE
        WRITE(*,*) "File ",'"',TRIM(filename),'"', " doesn't exist."
      END IF
      CLOSE(79)
! configuration displacement between different states
      IF(nstate==2) THEN
        ALLOCATE(dx(1:1,1:natom*3))
        DO i=1,natom
          dx(1,i*3-2)=xcoor(2,i*3-2)-xcoor(1,i*3-2)
          dx(1,i*3-1)=xcoor(2,i*3-1)-xcoor(1,i*3-1)
          dx(1,i*3)=xcoor(2,i*3)-xcoor(1,i*3)
        END DO
      ELSE IF(nstate==3) THEN
        ALLOCATE(dx(1:3,1:natom*3))   ! for more excited states!
      ELSE
        WRITE(*,*) "Check the "//TRIM(filename)//" file!"
      END IF
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
!!! calculate normal coordinate displacements: dqp !!!!!!!!!!!!!!!!!!!!!!
      dq=0.d0  ! unit: AMU^(1/2)*Angstrom
      dqp=0.d0 ! unit: Angstrom
      DO i=1,natom*3-6
        DO j=1,natom*3
          dqp(i)=dqp(i)+lmatrix(i,j)*dx(1,j)
        END DO
          dq(i)=dqp(i)*DSQRT(redmass(i))    ! unit: AMU^(1/2)*Angstrom
      END DO
!
! The following algorithm is used by C. K. Lin !!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculated mass weighted l matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     DO j=1,natom
!       lmw(1:natom*3-6,j*3-2)=lmatrix(1:natom*3-6,j*3-2)* &
!    &  DSQRT(n2m(0,mass(j)))
!       lmw(1:natom*3-6,j*3-1)=lmatrix(1:natom*3-6,j*3-1)* &
!    &  DSQRT(n2m(0,mass(j)))
!       lmw(1:natom*3-6,j*3)=lmatrix(1:natom*3-6,j*3)* &
!    &  DSQRT(n2m(0,mass(j)))
!     END DO
!!! re-normalize the lmw matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     f_norm=0.D0
!     f_norm(:)=SUM(lmw**2,DIM=2)
!     DO i=1,natom*3-6
!       DO j=1,natom*3
!         lmw(i,j)=lmw(i,j)/DSQRT(f_norm(i))
!       END DO
!     END DO
!!! calculate normal coordinate displacements: dq !!!!!!!!!!!!!!!!!!!!!!!
!     dq=0.d0
!     DO i=1,natom
!       dx(1,i*3-2)=dx(1,i*3-2)*DSQRT(n2m(0,mass(i)))
!       dx(1,i*3-1)=dx(1,i*3-1)*DSQRT(n2m(0,mass(i)))
!       dx(1,i*3)=dx(1,i*3)*DSQRT(n2m(0,mass(i)))
!     END DO
!! unit of dq: AMU^(1/2)*Angstrom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     DO i=1,natom*3-6
!       DO j=1,natom*3
!         dq(i)=dq(i)+lmw(i,j)*dx(1,j)
!       END DO
!     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      RETURN
      END SUBROUTINE Delta_Q
      END MODULE Delta_NMC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MODULE: FC_factor, Franck-Condon factor                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE FC_factor
      USE global
      USE Delta_NMC
      IMPLICIT NONE
!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Subroutine: FC_f
!    hrf(:): Huang-Rhys factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FC_f(natom,freq,dq,fscale,nu_max,fc)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom,nu_max
      REAL(KIND=8),INTENT(IN) :: fscale
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: freq,dq
      REAL(KIND=8),INTENT(OUT) :: fc(0:nu_max,1:natom*3-6)
! Local variables
      REAL(KIND=8),DIMENSION(natom*3-6) :: hrf,freq_au,dq_au
      INTEGER :: i,nu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the Displaced Osillator, the Huang-Rhys factor is computed as     !
! follows:                                                              !
!   S=(Beta*dqp^2)/2,                                                   !
!   Beta=(redmass*2Pi*Frequency)/hbar                                   !
! Note that unit of dqp(:) is Angstrom, and unit of freq(:) is cm-1, so !
!   Frequency=freq(:)*c,                                                !
! where c= is velocity of light. So,                                    !
!   S=(Pi*freq(:)*c*dq^2)/hbar                                          !
! The relation dq(:)=dqp(:)*redmass^(1/2) is used.                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO i=1,natom*3-6
        dq_au(i)=dq(i)*DSQRT(AMU2G/Me2G)/Bohr2A ! from AMU^(1/2)*Angstrom to (au)^(1/2)*au
        freq_au(i)=freq(i)*(clight*100.d0)*TAU  ! freq(:), from cm-1 to au(Time)
        hrf(i)=PI*freq_au(i)*(dq_au(i)**2)      ! hrf(i), unit in au, hbar=1.d0
! To calculate the Franck-Condon factor using the hrs(:)
        IF(hrf(i)<1.0d-25) hrf(i)=0.d0
        DO nu=0,nu_max
          fc(nu,i)=(hrf(i)**nu)*DEXP(-hrf(i))/factorials(nu)
          WRITE(*,*) nu,hrf(i),fc(nu,i)
        END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE FC_f
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Factorials: n!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION factorials(nu)
      IMPLICIT NONE
      INTEGER :: nutemp,ftemp
      INTEGER, INTENT(IN) :: nu
!
      ftemp=1
      nutemp=nu
      IF(nutemp<0) THEN
        WRITE(*,*) "Quantum number cannot be less than zero!"
      ELSE IF(nutemp<=1) THEN
        factorials=1
      ELSE
        DO WHILE(nutemp>1)
          ftemp=ftemp*nutemp
          nutemp=nutemp-1
        END DO
        factorials=ftemp
      END IF
      RETURN
      END FUNCTION factorials
!
      END MODULE FC_factor
