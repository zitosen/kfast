!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! gterm: calculate the G_l(t, T) terms                                  !
!                                                                       !
!             G_l=g_l \Prod_{j\=l}g_j                                   !
!                                                                       !
! References: (1) Hayashi, M., Lin, S.H., Raschke, M.B., Shen, Y.R.,    !
!                 2002. J. Phys. Chem. A 106, 2271–2282.                !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE spec_terms
      USE constants
!
      CONTAINS
      SUBROUTINE gterm(natom,freq,dq,sval,ntot,tstep,g_capital)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom,ntot
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: freq,dq,sval
      REAL(KIND=8),DIMENSION(natom*3-6,ntot+1),INTENT(OUT) :: g_capital
      REAL(KIND=8),DIMENSION(natom*3-6) :: wbolt
      COMPLEX(KIND=16),DIMENSION(natom*3-6,ntot+1) :: g_l,g_j
      REAL(KIND=8) :: tstep,t_kelvin
      INTEGER :: i,l,j
!
! (1) calculate Boltzmann weight: wbolt=(e^(hbar*omega/k_B*T)-1)^-1
      DO i=1,natom*3-6
        wbolt(i)=( EXP( h_Planck*freq(i)*clight*100 /     &
     &  (k_Boltzman*t_kelvin) ) -1 )**(-1)
      END DO
! (2) calculate the g_l term of the target mode l
!     note in complex exponential factor eiwt=e^(iwt), t is in fs, i.e.,
!     freq [cm-1] * clight [m/sec] * 100 ---> freq' [sec-1]
!     freq * clight * 1.0d-13 [fs-1]
      DO i=0,ntot
        DO l=1,natom*3-6
          g_l(l,i) = -0.5*dq(l)  &
             *(1-EXP(-cj*tstep*i *freq(l)*clight*1.0d-13)) &
     &       *EXP(-sval(l)*( (1+2*wbolt(l))                &
     &       -wbolt(l)*EXP(cj *tstep*i*freq(l)*clight*1.0d-13)  &
     &       -(1+wbolt(l))*EXP(-cj*tstep*i *freq(l)*clight*1.0d-13) ))
        END DO
      END DO
!
      END SUBROUTINE gterm
      END MODULE spec_terms
