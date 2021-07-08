!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! A code for doubly resonant sum-frequency generation (SFG)             !
! spectroscopy calculation.                                             !
!                                                                       !
!                                  Author: Z. Shen                      !
!                                  Email : shenzt@vip.henu.edu.cn       !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM dsfg_main
      USE constants
      USE nc_shift
      USE huang_rhys
      USE spec_terms 
      IMPLICIT NONE
      REAL(KIND=8) :: omega1,omega2,omega3,tstep,t_kelvin
      REAL(KIND=8),ALLOCATABLE :: freq(:),dq(:),sval(:)
      COMPLEX(KIND=8),ALLOCATABLE :: g_capital(:)
      INTEGER :: natom,ntot,l_mode,i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Initialization and Configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      natom=64                          ! Total atom number
      omega1=1.55956d0                  ! eV
      omega2=3.16834d0                  ! eV
      omega1=omega1/AU2ev               ! au
      omega2=omega2/AU2ev               ! au
      omega3=omega1+omega2              ! au
      ntot=5000                         ! total steps of time t
      tstep=1.d0                        ! time step, in femtosecond
      t_kelvin=300.d0                   ! Kelvin temperature
      l_mode=15                         ! the target mode
!
      ALLOCATE(freq(natom*3-6))
      ALLOCATE(dq(natom*3-6))
      ALLOCATE(sval(natom*3-6))
      ALLOCATE(g_capital(ntot+1))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now calculate Delta Normal Mode Coordinate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      CALL Delta_Q(natom,dq) 
!      WRITE(*,*) dq
!
! read the vibrational frequencies freq(:)
      OPEN(unit=77,file='freq.dat',status='old')
      DO i=1,natom*3-6
        READ(77,*) freq(i)    ! in cm-1
      END DO
      CLOSE(77)
!
! calculate Huang-Rhys factors
      CALL hrf(natom,dq,freq,sval)
!      WRITE(*,*) sval
!
! calculate G terms
!      WRITE(*,*) freq
      CALL gterm(natom,freq,dq,sval,ntot,l_mode,tstep,t_kelvin, &
     &                 g_capital)
 
      DO i=1,ntot+1
!        WRITE(*,FMT=10) REAL(g_capital(i)),' + i ', AIMAG(g_capital(i))
        WRITE(*,*) g_capital(i)
      ENDDO
 10   FORMAT(ES20.8, A, ES20.8)

!
      WRITE(*,'(A20)') "THANK GOD! ALL DONE!"
!
      END
