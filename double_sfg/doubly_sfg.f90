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
      USE NC_shift
      USE huang_rhys 
      IMPLICIT NONE
      REAL(KIND=8) :: omega1,omega2,omega3
      REAL(KIND=8),ALLOCATABLE :: freq(:),dq(:),sval(:),g_capital(:)
      INTEGER :: natom
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
!
      ALLOCATE(freq(natom*3-6))
      ALLOCATE(dq(natom*3-6))
      ALLOCATE(sval(natom*3-6))
      ALLOCATE(g_capital(natom*3-6))
!
! Now calculate Delta Normal Mode Coordinate
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
      CALL hrf(natom,dq,sval)
      WRITE(*,*) sval
!
! calculate G terms
      CALL gterm(natom,freq,dq,sval,g_capital)



!
      WRITE(*,'(A20)') "THANK GOD! ALL DONE!"
!
      END
