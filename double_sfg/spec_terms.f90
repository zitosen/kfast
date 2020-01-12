!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     gterm: calculate the G_l(t, T) terms                              !
!                                                                       !
!             G_l=g_l \Prod_{j\=l}g_j                                   !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE spec_terms
      USE constants
!
      CONTAINS
      SUBROUTINE gterm(natom,freq,sval,g_capital)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: dq
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(OUT) :: sval
      REAL(KIND=8),DIMENSION(natom*3-6) :: freq
      INTEGER :: i
!
! calculate 
      END SUBROUTINE gterm
      END MODULE spec_terms
