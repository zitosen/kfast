!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A simple molecular dynamics program
! refer to "D. Frenkel and B. Smit Understanding molecular simulation:
! from algorithms to applicaions, 2ed, Academic Press, 2009."
! Author: Zhitao Shen
! Email : shenzt@vip.henu.edu.cn
! Date  : 2021.1.25
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      PROGRAM md
      IMPLICIT NONE
      REAL(KIND=8) :: t,tmax,f,en,delta
!
      CALL init
      t=0
      DO WHILE (t.LT.tmax)
        CALL force(f,en)      ! determine the forces
        CALL integrate(f,en)  ! integrate equations of motion
        t=t+delta
        CALL sample           ! sample averages
      ENDDO
      STOP
      END
      
