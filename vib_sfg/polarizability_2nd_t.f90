!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate molecular vibrational second-order polarizability
! 2016.6.21, by zito
!
! Temperature effect is added
! 2017.10.21, by zito
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM POLARIZABILITY_2ND
      IMPLICIT NONE
      INTEGER, PARAMETER       :: dr=8,dc=16
      REAL(KIND=dr), PARAMETER :: amu2au=1.82288848426455D3,    &
                                  au2cm=2.1947463137054D5,      &
                                  au2eV=27.21138602D0,          &
                                  bohr2A=5.2917720859D-1,       &
                                  PI=3.141592653589793238462643D0, &
!                                  kb=1.38064852D-23 & ! J K^-1
                                  kb=8.6173303D-5  ! eV K^-1
! note here DQ should be consistent with 'normal_mode_shift.f90' code
      REAL(KIND=dr), PARAMETER :: DQ=1.0D-2
      INTEGER                  :: n_mode,ibegin,iend,nstep,i,j,k,l,m,n, &
                                  nu1,nu2,nu3,nu4,nu5,nu6,nu7,nu8,nu9,  &
                                  nu10,nu11,nu12,nu13,nu14,nu15,nu16,   &
                                  nu17,nu18,nu19,nu20,nu21,nu22,nu23,   &
                                  nu24,nu_max
      REAL(KIND=dr)            :: step,damp,photon2,angle_sum,angle_vis,&
                                  angle_ir,Lxx_vis,Lyy_vis,Lzz_vis,     &
                                  beta,temperature,temp_boltzmann
      COMPLEX(KIND=dc)         :: temp
      COMPLEX(KIND=dc), ALLOCATABLE, DIMENSION(:,:) &
                               :: polar_2nd,chi_eff
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:) &
                               :: freq,red_mass,Lxx_sum,Lyy_sum,Lzz_sum,&
                                  Lxx_ir,Lyy_ir,Lzz_ir
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) & 
                               :: polar_p,polar_m,dpolar,ddipole
      CHARACTER(LEN=100)       :: buffer,arg1
      CHARACTER(LEN=10)        :: subscript(27),sbuffer
!      
      CALL GETARG(1,arg1)
      arg1=TRIM(arg1)
      OPEN(UNIT=7,FILE=arg1,STATUS="OLD")
      READ(7,'(A)') buffer
      READ(7,*) n_mode,ibegin,iend,step,damp,temperature,nu_max
      WRITE(*,*) TRIM(buffer)
      WRITE(*,*) n_mode,ibegin,iend,step,damp
      ALLOCATE(freq(1:n_mode))
      ALLOCATE(red_mass(1:n_mode))
! frequency read in:
      READ(7,'(A)') buffer
      WRITE(*,*) TRIM(buffer)
      DO i=1,n_mode
        READ(7,*) freq(i),red_mass(i)
        WRITE(*,*) freq(i),red_mass(i)
      END DO
!
! dipole moment derivative read in:
      READ(7,'(A)') buffer
      WRITE(*,*) TRIM(buffer)
      READ(7,*) buffer
      ALLOCATE(ddipole(1:3,1:n_mode))
      DO i=1,n_mode
        READ(7,*) ddipole(1,i),ddipole(2,i),ddipole(3,i)
        WRITE(*,*) ddipole(1,i),ddipole(2,i),ddipole(3,i)
      END DO
! polarizability read in:
      READ(7,'(A)') buffer
      WRITE(*,*) TRIM(buffer)
      ALLOCATE(polar_p(1:9,1:n_mode))
      ALLOCATE(polar_m(1:9,1:n_mode))
      ALLOCATE(dpolar(1:9,1:n_mode))
! order of polarizability and its derivatives:
! column --> mode(1),mode(2),...; row --> xx,yy,zz,xy,xz,yz,yx,zx,zy
      DO i=1,n_mode
        READ(7,'(A)') buffer
        WRITE(*,*) TRIM(buffer)
        DO j=1,3
          READ(7,*) sbuffer,sbuffer,sbuffer,polar_p(j*3-2,i),polar_p(j*3-1,i),&
          polar_p(j*3,i)
          WRITE(*,*) polar_p(j*3-2,i),polar_p(j*3-1,i),polar_p(j*3,i)
        ENDDO
        READ(7,'(A)') buffer
        WRITE(*,*) TRIM(buffer)
        DO j=1,3
          READ(7,*) sbuffer,sbuffer,sbuffer,polar_m(j*3-2,i),polar_m(j*3-1,i),&
          polar_m(j*3,i)
          WRITE(*,*) polar_m(j*3-2,i),polar_m(j*3-1,i),polar_m(j*3,i)
        ENDDO
! now calculate polarizability derivative (angstrom^2*amu^(-1/2))
        dpolar(:,i)=(polar_p(:,i)-polar_m(:,i))/(2*DQ)
      ENDDO
      CLOSE(7)
!
! print out polarizability derivative
      OPEN(UNIT=7,FILE="ddpolar",STATUS="NEW")
      DO i=1,n_mode
        WRITE(7,*) "mode",i,"(angstrom^2*amu^(-1/2))"
        DO j=1,3
          WRITE(7,"(3F20.10)") dpolar(j*3-2,i),dpolar(j*3-1,i),dpolar(j*3,i)
        ENDDO
      ENDDO
      CLOSE(7)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now calculate second-order polarizability: polar_2nd
      nstep=INT((iend-ibegin+1)/step)
      ALLOCATE(polar_2nd(1:27,1:nstep))
      polar_2nd=(0.d0,0.d0)
      subscript=(/"xxx","xxy","xxz","yyx","yyy","yyz","zzx","zzy",      &
                  "zzz","xyx","xyy","xyz","xzx","xzy","xzz","yzx",      &
                  "yzy","yzz","yxx","yxy","yxz","zxx","zxy","zxz",      &
                  "zyx","zyy","zyz"/)
      freq(:)=freq(:)/au2cm ! tranfer 'freq' to 'au'
      m=0
      DO i=1,9    ! xx,yy,zz,xy,xz,yz,yx,zx,zy
        DO j=1,3  ! x,y,z
          m=m+1
          buffer="polar_2nd_"//TRIM(subscript(m))
          OPEN(UNIT=77,FILE=buffer,STATUS="NEW")
          WRITE(77,"(2A21)") "Frequencies(cm-1)    ","polar_2nd    "
          DO k=0,nstep-1
            photon2=ibegin+k*step
            temp=(0.d0,0.d0)
            DO l=1,n_mode
              IF (temperature.LT.1.0D-8) THEN
! for T=0 K, vibrational quantum number v=0, v'=1 (in a.u.) 
                temp=(0.d0,1.d0)*ddipole(j,l)*dpolar(i,l)               &
                *(bohr2A**(-2))*(amu2au**(-1/2))                        &
                *(2*(red_mass(l)*amu2au)*freq(l))**(-1)         &
                /((0.d0,1.d0)*(freq(l)*au2cm-photon2)/au2cm+damp/au2cm)
              ELSE
! for T>0 K, the temperature effect is considered by Boltzmann distribution
                beta=1.d0/(kb/au2eV*temperature) ! unit in au
                temp_boltzmann=1.d0
! mode=1
                IF(n_mode==1) THEN
                  DO nu1=0,nu_max
                    temp_boltzmann=DEXP(-nu1*freq(l)*beta) 
                    temp=temp+(0.d0,1.d0)*temp_boltzmann                &
                      *ddipole(j,l)*dpolar(i,l)                         &
                      *(bohr2A**(-2))*(amu2au**(-1/2))                  &
                      *((nu1+1)/(2*(red_mass(l)*amu2au)*freq(l)))       &
                      /((0.d0,1.d0)*(freq(l)*au2cm-photon2)/au2cm+damp/au2cm)
                  ENDDO
! mode=4
                ELSE IF(n_mode==4) THEN
                  DO nu1=0,nu_max
                  DO nu2=0,nu_max
                  DO nu3=0,nu_max
                  DO nu4=0,nu_max
                  IF(l==1) THEN 
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) 
                  ELSE IF(l==2) THEN
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta))
                  ELSE IF(l==3) THEN
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta))
                  ELSE IF(l==4) THEN
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)
                  ELSE
                    WRITE(*,*) "Check 'mode=4' case!"
                  ENDIF
                    temp=temp+(0.d0,1.d0)*temp_boltzmann                &
                      *ddipole(j,l)*dpolar(i,l)                         &
                      *(bohr2A**(-2))*(amu2au**(-1/2))                  &
                      *((nu1+1)/(2*(red_mass(l)*amu2au)*freq(l)))       &
                      /((0.d0,1.d0)*(freq(l)*au2cm-photon2)/au2cm+damp/au2cm)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
! mode=24
                ELSE IF(n_mode==24) THEN
                  DO nu1=0,nu_max
                  DO nu2=0,nu_max
                  DO nu3=0,nu_max
                  DO nu4=0,nu_max
                  DO nu5=0,nu_max
                  DO nu6=0,nu_max
                  DO nu7=0,nu_max
                  DO nu8=0,nu_max
                  DO nu9=0,nu_max
                  DO nu10=0,nu_max
                  DO nu11=0,nu_max
                  DO nu12=0,nu_max
                  DO nu13=0,nu_max
                  DO nu14=0,nu_max
                  DO nu15=0,nu_max
                  DO nu16=0,nu_max
                  DO nu17=0,nu_max
                  DO nu18=0,nu_max
                  DO nu19=0,nu_max
                  DO nu20=0,nu_max
                  DO nu21=0,nu_max
                  DO nu22=0,nu_max
                  DO nu23=0,nu_max
                  DO nu24=0,nu_max
                  SELECT CASE(l)
                  CASE(1)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta)) 
                  CASE(2)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(3)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(4)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(5)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(6)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(7)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(8)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(9)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(10)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(11)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(12)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(13)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(14)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(15)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(16)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(17)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(18)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(19)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(20)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(21)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(22)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(23)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta) &
                      *DEXP((-nu3*freq(24))*beta)*(1-DEXP(-freq(24)*beta))
                  CASE(24)
                    temp_boltzmann= &
                       DEXP((-nu1*freq(1))*beta)*(1-DEXP(-freq(1)*beta)) &
                      *DEXP((-nu2*freq(2))*beta)*(1-DEXP(-freq(2)*beta)) &
                      *DEXP((-nu3*freq(3))*beta)*(1-DEXP(-freq(3)*beta)) &
                      *DEXP((-nu4*freq(4))*beta)*(1-DEXP(-freq(4)*beta)) &
                      *DEXP((-nu2*freq(5))*beta)*(1-DEXP(-freq(5)*beta)) &
                      *DEXP((-nu3*freq(6))*beta)*(1-DEXP(-freq(6)*beta)) &
                      *DEXP((-nu4*freq(7))*beta)*(1-DEXP(-freq(7)*beta)) &
                      *DEXP((-nu2*freq(8))*beta)*(1-DEXP(-freq(8)*beta)) &
                      *DEXP((-nu3*freq(9))*beta)*(1-DEXP(-freq(9)*beta)) &
                      *DEXP((-nu4*freq(10))*beta)*(1-DEXP(-freq(10)*beta)) &
                      *DEXP((-nu2*freq(11))*beta)*(1-DEXP(-freq(11)*beta)) &
                      *DEXP((-nu3*freq(12))*beta)*(1-DEXP(-freq(12)*beta)) &
                      *DEXP((-nu4*freq(13))*beta)*(1-DEXP(-freq(13)*beta)) &
                      *DEXP((-nu2*freq(14))*beta)*(1-DEXP(-freq(14)*beta)) &
                      *DEXP((-nu3*freq(15))*beta)*(1-DEXP(-freq(15)*beta)) &
                      *DEXP((-nu4*freq(16))*beta)*(1-DEXP(-freq(16)*beta)) &
                      *DEXP((-nu2*freq(17))*beta)*(1-DEXP(-freq(17)*beta)) &
                      *DEXP((-nu3*freq(18))*beta)*(1-DEXP(-freq(18)*beta)) &
                      *DEXP((-nu4*freq(19))*beta)*(1-DEXP(-freq(19)*beta)) &
                      *DEXP((-nu2*freq(20))*beta)*(1-DEXP(-freq(20)*beta)) &
                      *DEXP((-nu3*freq(21))*beta)*(1-DEXP(-freq(21)*beta)) &
                      *DEXP((-nu4*freq(22))*beta)*(1-DEXP(-freq(22)*beta)) &
                      *DEXP((-nu2*freq(23))*beta)*(1-DEXP(-freq(23)*beta)) &
                      *DEXP((-nu3*freq(24))*beta)
                  CASE DEFAULT
                    WRITE(*,*) "Check 'mode=24'!"
                  END SELECT
                    temp=temp+(0.d0,1.d0)*temp_boltzmann                &
                      *ddipole(j,l)*dpolar(i,l)                         &
                      *(bohr2A**(-2))*(amu2au**(-1/2))                    &
                      *((nu1+1)/(2*(red_mass(l)*amu2au)*freq(l)))       &
                      /((0.d0,1.d0)*(freq(l)*au2cm-photon2)/au2cm+damp/au2cm)
                  ENDDO ! 1
                  ENDDO 
                  ENDDO
                  ENDDO
                  ENDDO ! 5
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO ! 10
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO ! 15
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO ! 20
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO ! 24
                ELSE
                  WRITE(*,*) "check 'mode selection condition!'"
                ENDIF ! condition for mode number
              END IF  ! condition for temperature
! sum of vibrational modes
              polar_2nd(m,k+1)=polar_2nd(m,k+1)+temp
!
! Note that the subscript order of second-order susceptibility in Lin's theory
! is different from other SFG theories, such as Y. Shen's formula.
! For example, polar_2nd_xyz is related to ddipole_y * dpolar_xz in Lin's 
! theory, while it is related to dpolar_xy * ddipole_z in Y. Shen's formula.
! This is because we specified different frequencies (omega1 or omega2) as IR light
! in the IR-visible SFG. Here we specify omega2 is IR light, and the subscript 
! order of second-order susceptibility is consistent with Y. Shen's formula.
! The order of polar_2nd (from m=1 to m=27) is:
! xxx xxy xxz yyx yyy yyz zzx zzy zzz
! xyx xyy xyz xzx xzy xzz yzx yzy yzz
! yxx yxy yxz zxx zxy zxz zyx zyy zyz
!
            END DO
! Here we print the polar_2nd
!            WRITE(77,"(F16.4,F20.12)") photon2,ABS(polar_2nd(m,k+1))
            WRITE(77,FMT=10) photon2,polar_2nd(m,k+1)
   10       FORMAT(F16.4,F20.12,' + ',F20.12)
          END DO
          CLOSE(77)
        END DO 
      END DO


!
      DEALLOCATE(freq)
      DEALLOCATE(red_mass)
      DEALLOCATE(ddipole)
      DEALLOCATE(polar_p)
      DEALLOCATE(polar_m)
      DEALLOCATE(dpolar)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now calculate the effective second-order susceptibilities for
! CH3OH/TiO2(110) system, which has C2v sysmetry. The lab frame x'y'z' 
! is obtained by rotating the crystallographic frame xyz with 225 degree in the plane.
! Refer to Liu et al., J. Phys. Chem. C 119 23486 (2015).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(UNIT=7,FILE="Fresnel_factor",STATUS="OLD")
      READ(7,'(A)') buffer
      READ(7,*) angle_sum,angle_vis,angle_ir
      angle_sum=angle_sum/180.D0*PI
      angle_vis=angle_vis/180.D0*PI
      angle_ir=angle_ir/180.D0*PI
      READ(7,'(A)') buffer
      READ(7,*) Lxx_vis,Lyy_vis,Lzz_vis
      READ(7,'(A)') buffer
      ALLOCATE(Lxx_sum(1:nstep))
      ALLOCATE(Lyy_sum(1:nstep))
      ALLOCATE(Lzz_sum(1:nstep))
      ALLOCATE(Lxx_ir(1:nstep))
      ALLOCATE(Lyy_ir(1:nstep))
      ALLOCATE(Lzz_ir(1:nstep))
      ALLOCATE(chi_eff(1:4,1:nstep))
      DO i=1,nstep
        READ(7,*) Lxx_sum(i),Lyy_sum(i),Lzz_sum(i),Lxx_ir(i),Lyy_ir(i), &
        Lzz_ir(i)
      END DO
      CLOSE(7)
!
      OPEN(UNIT=77,FILE="polar_eff",STATUS="NEW")
      WRITE(77,"(A15,A30)") "Freq(cm-1)  ","ssp,    sps,   pss,   ppp"
! Here chi_eff(1,:), chi_eff(2,:),..., correspond to ssp, sps, pss, ppp.
      chi_eff(1,:)=0.5d0*Lyy_sum(:)*Lyy_vis*Lzz_ir(:)*DSIN(angle_ir)  &
                   *(polar_2nd(3,:)+polar_2nd(6,:))
      chi_eff(2,:)=0.5d0*Lyy_sum(:)*Lzz_vis*Lyy_ir(:)*DSIN(angle_vis) &
                   *(polar_2nd(13,:)+polar_2nd(17,:))
      chi_eff(3,:)=0.5d0*Lzz_sum(:)*Lyy_vis*Lyy_ir(:)*DSIN(angle_sum) &
                   *(polar_2nd(22,:)+polar_2nd(26,:))
      chi_eff(4,:)=-0.5d0*Lxx_sum(:)*Lxx_vis*Lzz_ir(:)*DCOS(angle_sum) &
                   *DCOS(angle_vis)*DSIN(angle_ir) &
                   *(polar_2nd(3,:)+polar_2nd(6,:))      &
                   -0.5d0*Lxx_sum(:)*Lzz_vis*Lxx_ir(:)*DCOS(angle_sum) &
                   *DSIN(angle_vis)*DCOS(angle_ir) &
                   *(polar_2nd(13,:)+polar_2nd(17,:))    &
                   +0.5d0*Lzz_sum(:)*Lxx_vis*Lxx_ir(:)*DSIN(angle_sum) &
                   *DCOS(angle_vis)*DCOS(angle_ir) &
                   *(polar_2nd(22,:)+polar_2nd(26,:))    &
                   +Lzz_sum(:)*Lzz_vis*Lzz_ir(:)*DSIN(angle_sum) &
                   *DSIN(angle_vis)*DSIN(angle_ir) &
                   *polar_2nd(9,:)
      DO i=1,nstep
      WRITE(77,'(4F20.12)') ABS(chi_eff(1,i)),ABS(chi_eff(2,i)), &
                       ABS(chi_eff(3,i)),ABS(chi_eff(4,i))
      END DO
      CLOSE(77)
      OPEN(77,FILE="chi_ssp",STATUS="NEW")
        DO i=1,nstep
        WRITE(77,'(2F20.12)') chi_eff(1,i)
        END DO
      CLOSE(77)
      OPEN(77,FILE="chi_sps",STATUS="NEW")
        DO i=1,nstep
        WRITE(77,'(2F20.12)') chi_eff(2,i)
        END DO
      CLOSE(77)
      OPEN(77,FILE="chi_pss",STATUS="NEW")
        DO i=1,nstep
        WRITE(77,'(2F20.12)') chi_eff(3,i)
        END DO
      CLOSE(77)
      OPEN(77,FILE="chi_ppp",STATUS="NEW")
        DO i=1,nstep
        WRITE(77,'(2F20.12)') chi_eff(4,i)
        END DO
      CLOSE(77)
!
      DEALLOCATE(polar_2nd)
      DEALLOCATE(Lxx_sum)
      DEALLOCATE(Lyy_sum)
      DEALLOCATE(Lzz_sum)
      DEALLOCATE(Lxx_ir)
      DEALLOCATE(Lyy_ir)
      DEALLOCATE(Lzz_ir)
      DEALLOCATE(chi_eff)     
      END PROGRAM

