!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Sum-frequency generation (SFG) spectroscopy calculation program.      !
!                                                                       !
!                    2015/12/15  by Z. Shen                             !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM sfg_main
      USE ir_raman_tensor
      USE FC_factor
      IMPLICIT NONE
      REAL(KIND=8) :: alpha(1:3,1:3),fscale
      REAL(KIND=8),ALLOCATABLE :: freq(:),dq(:),fc(:,:)
      INTEGER :: natom,nstate_max,nstate,nu_max
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Initialization and Configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      WRITE(6,*) " The resonant state is (like 1,2,...)"
!      READ(5,*) i_resonant
!      WRITE(6,*) " The number of excited states considered:"
!      READ(5,*) nstate
!      WRITE(6,*) " The omega_1 and omega_2 (omega_1<omeag_2, eV):"
!      READ(5,*) omega(1),omega(2)
!

      natom=10                          ! Total atom number
      nstate_max=200
      i_resonant=1                      ! number of resonant state
      omega(1)=1.55956d0                ! eV
      omega(2)=3.16834d0                ! eV
      omega(1)=omega(1)/au2ev           ! au
      omega(2)=omega(2)/au2ev           ! au
      omega(3)=omega(1)+omega(2)        ! au
! files incluing data
      file_TDM="S0Sn_TEDM_cis_6-31gd.log"
      file_E="S0Sn_TEDM_cis-d_6-31gd.log"
!
! Now calculate polarizability-like terms, alpha !!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(unit=11,file='alpha.dat',status='unknown')
      WRITE(11,'(A32)') " nstate   alpha_xy    alpha_yx"
      DO nstate=0,nstate_max
        CALL tensor_cal(alpha,nstate)
        WRITE(11,'(1X,I4,2F18.10)') nstate,alpha(1,2),alpha(2,1)
      END DO
      CLOSE(11)
!
! Now calculate Delta Normal Mode Coordinate !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(freq(natom*3-6))
      ALLOCATE(dq(natom*3-6))
      CALL Delta_Q(natom,freq,dq) 
!
!
! Now calculate Franck-Condon factor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      fscale=0.96d0    ! The caled factor to harmonic vibration frequencies
      nu_max=6         ! The highest value of vibrational quantum number for
!                      !  excited electronic state.
      ALLOCATE(fc(0:nu_max,1:natom*3-6))
      CALL FC_f(natom,freq,dq,fscale,nu_max,fc)
!      WRITE(*,*) fc(:,2)
!
      WRITE(*,'(A20)') "THANK GOD! ALL DONE!"
!
      END
