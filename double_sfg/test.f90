      PROGRAM test
      IMPLICIT NONE
      INTEGER, PARAMETER :: natom=64, ntot=5000 
      COMPLEX(KIND=16),DIMENSION(natom*3-6,ntot+1) :: g_low
!
      g_low(:,:)=(1.d0,0.d0)
      WRITE(*,*) g_low
      END
