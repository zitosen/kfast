! Assuming a simple cubic lattice 
      SUBROUTINE init
      IMPLICIT NONE
      REAL(KIND=8) :: sumv,sumv2
      REAL(KIND=8), ALLOCATABLE :: x(:),lattice_pos(:)
! npart, number of particles
      INTEGER :: npart,i
!
      sumv=0
      sumv2=0
      DO i=1,npart
        x(i)=lattice_pos(i) ! place the particles on a lattice
        v(i)=(RANF()-0.5)   ! give random velocities
