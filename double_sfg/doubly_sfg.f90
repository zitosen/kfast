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
      IMPLICIT NONE
      REAL(KIND=8) :: omega1,omega2,omega3
      REAL(KIND=8),ALLOCATABLE :: freq(:),dq(:)
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
! Now calculate Delta Normal Mode Coordinate !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      ALLOCATE(freq(natom*3-6))
      ALLOCATE(dq(natom*3-6))
      CALL Delta_Q(natom,dq) 
      WRITE(*,*) dq

!
      WRITE(*,'(A20)') "THANK GOD! ALL DONE!"
!
      END
