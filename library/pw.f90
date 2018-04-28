!-------------------------------------------------------------------
      program plane_wave
!-------------------------------------------------------------------
! Purpose: Solve one-dimensional Schrodinger eqation H|psi> = E|psi>,
!          where H=p^2/2m+V(x), V(x) is a finite (symmetric) potential 
!          well, i.e., 
!          V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2,
!          |psi> is represented by a plane-wave basis set
!          b_n=1/sqrt(a)*e^(i k_n x)
! Units  : atomic Rydberg units, hbar^2/2m = 1
! Library: lapack, dsyev

      implicit none
      integer, parameter :: dp = selected_real_kind(14,200)
      real(dp), parameter :: pi=3.14159265358979_dp
      integer :: n, npw
      real(dp) :: v0, a, b
      real(dp), allocatable :: kn(:), e(:), h(:,:), work (:)
      real(dp) :: x, dx, norm, prob
      complex(dp) :: f
      integer :: i, j, nr, lwork, info
!
! Input data
!
      write (*,"('Parameters for potential well: V_0, b > ',$)")
      read (*,*) v0, b
      if (v0 <= 0.0_dp .or. b <= 0.0_dp) stop ' wrong input parameters '
      write (*,"('   V_0, b =',2f10.4)") v0, b
!
! Plane waves between -a/2 < x < a/2, k_i= +- 2*pi*i/a, i=0,1,...,n
!
      write (*,"('Parameters for plane waves: a, n > ',$)")
      read (*,*) a, n
      if ( a <= b .or. n <= 0) stop ' wrong input parameters '
      write (*,"('a, n=',f8.4,i6)") a, n
      npw = 2*n+1
      allocate (kn(npw), e(npw), work(3*npw), h(npw,npw) )
!
! Assign values of k_n: n=0,+1,-1,+2,-2, etc
!
      kn(1) = 0.0_dp
      do i=2,npw-1,2
        kn(i)   = (i/2)*2.0_dp*pi/a
        kn(i+1) =-(i/2)*2.0_dp*pi/a
      end do
!
! Assign values of the matrix elements of the hamiltonian on the
! plane-wave basis,
!   H_ij = P_ij + V_ij
!   P_ij = (hbar * k_i)^2 / (2m) * delta_ij
!
!          -V_0 sin[b(k_i-k_j)/2]
!   V_ij =  --- -----------------, i /= j
!            a    (k_i-k_j)/2
!
!   V_ii = -(V_0 * b) / a
!
      h(:,:) = 0.0_dp
      do i=1,npw
        do j=1,npw
          if ( i ==j ) then
            h(i,j) = kn(i)**2 - v0*b/a
          else
            h(i,j) = -v0/a * sin( (kn(j)-kn(i))*b/2.0_dp ) /  &
           (kn(j)-kn(i))*2.0_dp
          end if
        end do
      end do
!
! Solution [expansion coefficients are stored into h(j,i)
! j=basis function index, i= eigenvalue index]
!
      lwork = 3*npw
! The "leading dimension of array" h is its first dimension
! For dynamically allocated matrices, see the "allocate" command
      call dsyev ( 'V', 'U', npw, h, npw, e, work, lwork, info )
      if (info /= 0) stop 'H-matrix diagonalization failed '
!
      write (*,"('   Lowest eigenvalues:',3f12.6)") (e(i),i=1,3)
! output the ground state wavefunction
      open (7,file='gs_wfc.out',status='unknown',form='formatted')
      dx = 0.01_dp
      nr = nint(a/2.0_dp/dx)
      norm = 0.0_dp
      do i=-nr, nr
        x = dx*i
        f = 0.0_dp
        do j=1,npw
          f = f + h(j,1)*exp((0.0,1.0)*kn(j)*x)/sqrt(a)
        end do
        prob = f*conjg(f)
        norm = norm + prob*dx
        write(7,'(f12.6,3f10.6)') x, prob, f
      end do 
! verify normalization:
      write (*,"('   norm: ',f12.6)") norm
      close(7)
      deallocate ( h, work, e, kn)
      end program plane_wave 
