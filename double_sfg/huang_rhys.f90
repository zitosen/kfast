!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     huang_rhys: calculate the Huang-Rhys factors S_i by               !
!                                                                       !
!            S_i = omega_i*dq_i^2/(2*hbar)                              !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE huang_rhys
      USE constants
!
      CONTAINS
      SUBROUTINE hrf(natom,dq,freq,sval)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: natom
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: dq
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(IN) :: freq
      REAL(KIND=8),DIMENSION(natom*3-6),INTENT(OUT) :: sval
      INTEGER :: i
!
! calculate Huang-Rhys factors
! freq [cm-1] * clight [cm/sec] ---> freq' [sec-1]
!      freq * clight * 100
! dq [AMU^(1/2)*Angstrom] * AMU2kg^(1/2) * 10^(-10)  ---> kg^(1/2)*m
!      dq**2 * AMU2kg * 1.d-20
! hbar in J*sec
!      hbar=h_Planck/(2*PI)
!       
      DO i=1,natom*3-6
        sval(i)=(2*PI**2*freq(i)*clight*100     &          
     &           * dq(i)**2 * AMU2kg * 1.d-20)  & 
     &           * h_Planck**(-1)
      END DO
!
      END SUBROUTINE hrf
      END MODULE huang_rhys
