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
! l_mode is the target mode
! ntot is the total time steps
      SUBROUTINE gterm(natom,freq,dq,sval,ntot,l_mode,tstep,t_kelvin, &
     &                 g_capital,omega1,omega2,d_capital)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom,ntot,l_mode
      REAL(KIND=8),INTENT(IN) :: tstep,t_kelvin,omega1,omega2
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: freq,dq,sval
      REAL(KIND=8),DIMENSION(natom*3-6) :: wbolt
      COMPLEX(KIND=8),DIMENSION(ntot+1),INTENT(OUT) :: g_capital
      COMPLEX(KIND=8),INTENT(OUT) :: d_capital
      COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: g_low  ! DIMENSION(ntot+1, natom*3-6) statement
!                                                             ! results in "Segmentation fault"
      INTEGER :: i,l
!
! (1) calculate Boltzmann weight: wbolt=(e^(hbar*omega/k_B*T)-1)^-1
!     h_Planck  : J*second
!     freq      : cm-1
!     k_Boltzman: J/degree
!     T         : Kelvin
      DO i=1,natom*3-6
        wbolt(i)=( EXP( (h_Planck/(2*PI))*freq(i)*clight*100 /     &
     &  (k_Boltzman*t_kelvin) ) -1 )**(-1)
!        WRITE(*,*) wbolt(i)
      END DO
! (2) calculate the g_l and g_j terms
!     note in complex exponential factor eiwt=e^(iwt), t is in fs, i.e.,
!     freq [cm-1] * clight [m/sec] * 100 ---> freq' [sec-1]
!     freq * clight * 1.0d-13 [fs-1]
      ALLOCATE(g_low(1:ntot+1, 1:natom*3-6))
      DO i=1,ntot+1
        DO l=1,natom*3-6
          IF(l==l_mode) THEN
! dq in unit AMU^(1/2)*Angstrom
            g_low(i,l) = (-0.5*dq(l))                                 &
     &       *(1-EXP(-cj*tstep*i *freq(l)*clight*1.0d-13))            &
     &       *EXP(-sval(l)*( (1+2*wbolt(l))                           &
     &       -wbolt(l)*EXP(cj *tstep*i*freq(l)*clight*1.0d-13)        &
     &       -(1+wbolt(l))*EXP(-cj*tstep*i *freq(l)*clight*1.0d-13) ))
          ELSE
            g_low(i,l) = EXP( -sval(l)*( (1+2*wbolt(l))               &
     &       -wbolt(l)*EXP(cj*tstep*i*freq(l)*clight*1.0d-13)         &
     &       -(1+wbolt(l))*EXP(-cj*tstep*i*freq(l)*clight*1.0d-13) ))
          END IF
        END DO
      END DO
! (3) calculate the G_l(t, T) terms
      g_capital=(1.d0,0.d0)
      DO i=1,natom*3-6
        g_capital(:)=g_capital(:)*g_low(:,i)
      END DO
! (4) calculate D_l(omega1,omega2,T)
      d_capital=(0.d0,0.d0)
      DO i=1,ntot+1,2
        EXP(-tstep*i(cj*() ))
      ENDDO





      DEALLOCATE(g_low)
      RETURN
      END SUBROUTINE gterm
      END MODULE spec_terms
