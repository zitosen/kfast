! The code to estimate tolerance factors of organic-inorganic hybrid
! perovskites.
! 2020.6.30, by zito (shenzt@vip.henu.edu.cn)
! 
      PROGRAM Tolerance_F
      IMPLICIT NONE
      INTEGER, PARAMETER       :: sr=4,sc=8,dr=8,dc=16
      REAL(KIND=dr), PARAMETER :: rho_step=0.01d0
      REAL(KIND=dr), PARAMETER :: r_Cs=0.167d0,r_FA=0.253d0,r_MA=0.217d0,&
                                  r_Pb=0.119d0,r_I=0.220d0,r_Br=0.196d0
      REAL(KIND=dr) :: x,y
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) :: toler_f
      INTEGER                  :: i,j,tot_rho
      
!
! Perovskite stoichiometric ratios: "Cs_x FA_(1-x) Pb I_3(1-y) Br_3y"
      OPEN(UNIT=7,FILE="tolerance.dat",STATUS="UNKNOWN")
      tot_rho=1.0/rho_step+1
      ALLOCATE(toler_f(1:tot_rho,1:tot_rho))
      DO i=1,tot_rho   ! for x
        DO j=1,tot_rho ! for y
          x=(i-1)*rho_step
          y=(j-1)*rho_step
          toler_f(i,j)=( (r_Cs*x + r_FA*(1-x)) + r_I*(1-y)              &
     &                                         + r_Br*y  )              &
     &                / ( SQRT(2.0)* ( r_Pb + r_I*(1-y)                 &
     &                                      + r_Br*y ) )
        WRITE(7,"(3F10.4)") x, y, toler_f(i,j)
        ENDDO
        WRITE(7,*) " "
      ENDDO
      CLOSE(7)
      END PROGRAM
