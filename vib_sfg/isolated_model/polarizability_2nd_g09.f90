! Calculate molecular vibrational second-order polarizability
! dipole and polarizability derivatives are obtained from G09
! 2017.5.15, by zito
!
      PROGRAM POLARIZABILITY_2ND_G09
      IMPLICIT NONE
      INTEGER, PARAMETER       :: sr=4,sc=8,dr=8,dc=16
      REAL(KIND=sr), PARAMETER :: num_den=1.0  ! molecular number density
      REAL(KIND=dr), PARAMETER :: amu2au=1.82288848426455D3,    &
     & au2cm=2.1947463137054D5, bohr2A=5.2917720859D-1, &
     & PI=3.14159265358979323846264338328d0
      INTEGER                  :: n_mode,ibegin,iend,nstep,i,j,k,l,m
      REAL(KIND=dr)            :: step,damp,photon1,photon2,beta,beta1,&
     & beta2,theta,psi
      COMPLEX(KIND=dc)         :: temp
      COMPLEX(KIND=dc), ALLOCATABLE, DIMENSION(:) &
                               :: chi_ssp,chi_sps,chi_pss,chi_ppp,&
     & chi_spp,chi_pps,chi_psp
      COMPLEX(KIND=dc), ALLOCATABLE, DIMENSION(:,:) &
                               :: polar_2nd,polar_2nd_ave
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:) &
                               :: freq,red_mass
      REAL(KIND=dr),DIMENSION(3) &
                               :: fresnel
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) & 
                               :: dpolar,ddipole
      CHARACTER(LEN=100)       :: buffer,arg1
      CHARACTER(LEN=10)        :: subscript(27),sbuffer
!     
      WRITE(*,*) "Number of normal modes considered: "
      READ(*,*) n_mode 
      ALLOCATE(freq(1:n_mode))
      ALLOCATE(red_mass(1:n_mode))
! frequency read in:
      OPEN(UNIT=7,FILE="freq.dat",STATUS="OLD")
      DO i=1,n_mode
        READ(7,*) freq(i)
!        WRITE(*,*) "Read freqency, Done!"
      END DO
      CLOSE(7)
! reduced mass read in:
      OPEN(UNIT=7,FILE="reducedmass.dat",STATUS="OLD")
      DO i=1,n_mode
        READ(7,*) red_mass(i)
!        WRITE(*,*) "Read reduced mass, Done!"
      END DO
      CLOSE(7) 
! dipole moment derivatives read in:
      ALLOCATE(ddipole(1:3,1:n_mode))
      OPEN(UNIT=7,FILE="ddipole.dat",STATUS="OLD")
      DO i=1,n_mode
        READ(7,*) ddipole(1,i),ddipole(2,i),ddipole(3,i)
!        WRITE(*,*) "Read ddipole, Done!"
      END DO
      CLOSE(7)
! polarizability derivatives read in:
!   1-9 for  xx xy xz  yx yy yz  zx zy zz
      ALLOCATE(dpolar(1:9,1:n_mode))
      OPEN(UNIT=7,FILE="dpolar.dat",STATUS="OLD")
      DO i=1,n_mode
        DO j=1,3
          READ(7,*) dpolar(1+3*(j-1),i),dpolar(2+3*(j-1),i),dpolar(3+3*(j-1),i)
!          WRITE(*,*) "Read dpolar, Done!"
        ENDDO
      ENDDO
      CLOSE(7)
!
! calculate second-order polarizability: polar_2nd
      WRITE(*,*) "Now calculate second-order polarizability:"
      WRITE(*,*) "Freqency range and step of second-order polarizability (cm-1), &
     &            eg, 2000 4000 1 "
      READ(*,*) ibegin,iend,step
      WRITE(*,*) "Damping constant Gamma (cm-1):"
      READ(*,*) damp
      nstep=INT((iend-ibegin+1)/step)
      ALLOCATE(polar_2nd(1:27,1:nstep))
      polar_2nd=(0.d0,0.d0)
      subscript=(/"xxx","xxy","xxz","xyx","xyy","xyz","xzx","xzy",&
                  "xzz","yxx","yxy","yxz","yyx","yyy","yyz","yzx",&
                  "yzy","yzz","zxx","zxy","zxz","zyx","zyy","zyz",&
                  "zzx","zzy","zzz"/)
      m=0
      DO i=1,9    ! xx xy xz  yx yy yz  zx zy zz
        DO j=1,3  ! x,y,z
          m=m+1
          buffer="polar_2nd_"//TRIM(subscript(m))
          OPEN(UNIT=77,FILE=buffer,STATUS="UNKNOWN")
          WRITE(77,"(2A21)") "Frequencies(cm-1)    ",buffer
          DO k=0,nstep-1
            photon2=ibegin+k*step
            temp=(0.d0,0.d0)
            DO l=1,n_mode
              ! vibrational quantum number v=0, v'=1; T=0 K
              ! convert to a.u. 
              temp=(0.d0,1.d0)*ddipole(j,l)*dpolar(i,l)                 &
              *(bohr2A**(-4))*(amu2au**(-1))                            &
              *(2*(red_mass(l)*amu2au)*(freq(l)/au2cm))**(-1)           &
              /((0.d0,1.d0)*(freq(l)-photon2)/au2cm+damp/au2cm)
              ! sum of vibrational modes
              polar_2nd(m,k+1)=polar_2nd(m,k+1)+temp
            END DO
            WRITE(77,"(F10.2,E20.10)") photon2,ABS(polar_2nd(m,k+1))
          END DO
          CLOSE(77)
        END DO 
      END DO
      DEALLOCATE(freq)
      DEALLOCATE(red_mass)
      DEALLOCATE(ddipole)
      DEALLOCATE(dpolar)
! Orientional Averages through Euler Rotation Matrix
! Here assuming a random distribution in the in-plane rotation angle phi,
! i.e., yielding C_infinity macroscopic symmetry. There are seven unique elements
! of chi^(2), chi_zzz,  chi_zxx=chi_zyy, chi_xzx=chi_yzy, chi_xxz=chi_yyz,
! chi_xyz=-chi_yxz, chi_xzy=-chi_yzx, chi_zxy=-chi_zyx. Refer to A. J. Moad and 
! G. J. Simpson, J. Phys. Chem. B 2004, 108, 3548. (Note the Euler rotation matrix
! is in Z-y'-z convention, (X Y Z)^T=R*(x y z)^T)
! For C1 molecule, assuming delta function distributions of Euler angles theta
! and psi
      ALLOCATE(polar_2nd_ave(1:7,1:nstep))
      WRITE(*,*) "Euler angles theta (tilt) and psi (twist) in degree: "
      READ(*,*) theta,psi
      theta=theta/180.d0*PI
      psi=psi/180.d0*PI
! chi_ZZZ
      polar_2nd_ave(1,:)=num_den & 
     &*DCOS(theta)**3*polar_2nd(27,:)                        & !zzz
     &+DSIN(theta)*DSIN(psi)                                 &
     & *(polar_2nd(18,:)+polar_2nd(24,:)+polar_2nd(26,:))    & !yzz zyz zzy
     &-DSIN(theta)*DCOS(psi)                                 &
     & *(polar_2nd(9,:)+polar_2nd(21,:)+polar_2nd(25,:))     & !xzz zxz zzx
     &+DSIN(theta)**2*DCOS(theta)*DSIN(psi)**2               &
     & *(polar_2nd(15,:)+polar_2nd(17,:)+polar_2nd(23,:))    & !yyz yzy zyy
     &+DSIN(theta)**2*DCOS(theta)*DCOS(psi)**2               &
     & *(polar_2nd(3,:)+polar_2nd(7,:)+polar_2nd(19,:))      & !xxz xzx zxx
     &-DSIN(theta)**2*DCOS(theta)*DSIN(psi)*DCOS(psi)        &
     & *(polar_2nd(6,:)+polar_2nd(8,:)+polar_2nd(12,:)       & ! xyz xzy yxz
     & +polar_2nd(16,:)+polar_2nd(20,:)+polar_2nd(22,:))     & ! yzx zxy zyx
     &+DSIN(theta)**3*DSIN(psi)                              &
     & *(polar_2nd(2,:)+polar_2nd(4,:)+polar_2nd(10,:)       & ! xxy xyx yxx
     & -polar_2nd(18,:)-polar_2nd(24,:)-polar_2nd(26,:))     & ! yzz zyz zzy
     &+DSIN(theta)**3*DCOS(psi)                              &
     & *(polar_2nd(9,:)+polar_2nd(7,:)+polar_2nd(25,:)       & ! xzz zxz zzx
     & -polar_2nd(5,:)-polar_2nd(11,:)-polar_2nd(13,:))      & ! xyy yxy yyx
     &+DSIN(theta)**3*DSIN(psi)**3                           &
     & *(polar_2nd(14,:)-polar_2nd(2,:)-polar_2nd(4,:)       & ! yyy xxy xyx
     & -polar_2nd(10,:))                                      & ! yxx
     &+DSIN(theta)**3*DCOS(psi)**3                           &
     & *(-polar_2nd(1,:)+polar_2nd(5,:)+polar_2nd(11,:)      & ! xxx xyy yxy
     & +polar_2nd(13,:))                                        ! yyx
!
!      DO i=1,nstep
!        WRITE(*,*) polar_2nd_ave(1,i)
!      ENDDO
! chi_ZXX
      polar_2nd_ave(2,:)=num_den*0.5                         &
     &*DSIN(theta)**2*DCOS(theta)*polar_2nd(27,:)            & ! zzz
     &+DCOS(theta)*(polar_2nd(19,:)+polar_2nd(23,:))         & ! zxx zyy
     &-DSIN(theta)**2*DCOS(theta)*DSIN(psi)**2               &
     & *(polar_2nd(15,:)+polar_2nd(17,:)+polar_2nd(23,:))    & ! yyz yzy zyy
     &-DSIN(theta)**2*DCOS(theta)*DCOS(psi)**2               &
     & *(polar_2nd(3,:)+polar_2nd(7,:)+polar_2nd(19,:))      & ! xxz xzx zxx
     &+DSIN(theta)**2*DCOS(theta)*DSIN(psi)*DCOS(psi)        &
     & *(polar_2nd(6,:)+polar_2nd(8,:)+polar_2nd(12,:)       & ! xyz xzy yxz
     & +polar_2nd(16,:)+polar_2nd(20,:)+polar_2nd(22,:))     & ! yzx zxy zyx
     &+DSIN(theta)*DSIN(psi)                                 &
     & *(polar_2nd(14,:)+polar_2nd(10,:)-polar_2nd(24,:)     & ! yyy yxx zyz
     & -polar_2nd(26,:))                                     & ! zzy
     &+DSIN(theta)*DCOS(psi)                                 &
     & *(-polar_2nd(1,:)-polar_2nd(5,:)+polar_2nd(21,:)      & ! xxx xyy zxz
     & +polar_2nd(25,:))                                     & ! zzx
     &+DSIN(theta)**3*DSIN(psi)                              &
     & *(-polar_2nd(2,:)-polar_2nd(4,:)-polar_2nd(10,:)      & ! xxy xyx yxx
     & +polar_2nd(18,:)+polar_2nd(24,:)+polar_2nd(26,:))     & ! yzz zyz zzy
     &+DSIN(theta)**3*DCOS(psi)                              &
     & *(polar_2nd(5,:)+polar_2nd(11,:)+polar_2nd(13,:)      & ! xyy yxy yyx
     & -polar_2nd(9,:)-polar_2nd(21,:)-polar_2nd(25,:))      & ! xzz zxz zzx
     &+DSIN(theta)**3*DSIN(psi)**3                           &
     & *(-polar_2nd(14,:)+polar_2nd(2,:)+polar_2nd(4,:)      & ! yyy xxy xyx
     & +polar_2nd(10,:))                                     & ! yxx
     &+DSIN(theta)**3*DCOS(psi)**3                           &
     & *(polar_2nd(1,:)-polar_2nd(5,:)-polar_2nd(11,:)       & ! xxx xyy yxy
     & -polar_2nd(13,:))                                       ! yyx
! chi_XZX
      polar_2nd_ave(3,:)=num_den*0.5                         &
     &*DSIN(theta)**2*DCOS(theta)*polar_2nd(27,:)            & ! zzz
     &+DCOS(theta)*(polar_2nd(7,:)+polar_2nd(17,:))          & ! xzx yzy
     &-DSIN(theta)**2*DCOS(theta)*DSIN(psi)**2               & 
     & *(polar_2nd(15,:)+polar_2nd(17,:)+polar_2nd(23,:))    & ! yyz yzy zyy
     &-DSIN(theta)**2*DCOS(theta)*DCOS(psi)**2               &
     & *(polar_2nd(3,:)+polar_2nd(7,:)+polar_2nd(19,:))      & ! xxz xzx zxx
     &+DSIN(theta)**2*DCOS(theta)*DSIN(psi)*DCOS(psi)        &
     & *(polar_2nd(6,:)+polar_2nd(8,:)+polar_2nd(12,:)       & ! xyz xzy yxz
     & +polar_2nd(16,:)+polar_2nd(20,:)+polar_2nd(22,:))     & ! yzx zxy zyx
     &+DSIN(theta)*DSIN(psi)                                 &
     & *(polar_2nd(14,:)+polar_2nd(4,:)-polar_2nd(18,:)      & ! yyy xyx yzz
     & -polar_2nd(26,:))                                     & ! zzy
     &+DSIN(theta)*DCOS(psi)                                 &
     & *(-polar_2nd(1,:)-polar_2nd(11,:)+polar_2nd(9,:)      & ! xxx yxy xzz
     & +polar_2nd(25,:))                                     & ! zzx
     &+DSIN(theta)**3*DSIN(psi)                              &
     & *(-polar_2nd(2,:)-polar_2nd(4,:)-polar_2nd(10,:)      & ! xxy xyx yxx
     & +polar_2nd(18,:)+polar_2nd(24,:)+polar_2nd(26,:))     & ! yzz zyz zzy
     &+DSIN(theta)**3*DCOS(psi)                              &
     & *(polar_2nd(5,:)+polar_2nd(11,:)+polar_2nd(13,:)      & ! xyy yxy yyx
     & -polar_2nd(9,:)-polar_2nd(21,:)-polar_2nd(25,:))      & ! xzz zxz zzx
     &+DSIN(theta)**3*DSIN(psi)**3                           &
     & *(-polar_2nd(14,:)+polar_2nd(2,:)+polar_2nd(4,:)      & ! yyy xxy xyx
     & +polar_2nd(10,:))                                     & ! yxx
     &+DSIN(theta)**3*DCOS(psi)**3                           &
     & *(polar_2nd(1,:)-polar_2nd(5,:)-polar_2nd(11,:)       & ! xxx xyy yxy
     & -polar_2nd(13,:))                                       ! yyx
! chi_XXZ
      polar_2nd_ave(4,:)=num_den*0.5                         &
     &*DSIN(theta)**2*DCOS(theta)*polar_2nd(27,:)            & ! zzz
     &+DCOS(theta)*(polar_2nd(3,:)+polar_2nd(15,:)) & ! xxz yyz
     &-DSIN(theta)**2*DCOS(theta)*DSIN(psi)**2 & 
     & *(polar_2nd(15,:)+polar_2nd(17,:)+polar_2nd(23,:)) & ! yyz yzy zyy
     &-DSIN(theta)**2*DCOS(theta)*DCOS(psi)**2 &
     & *(polar_2nd(3,:)+polar_2nd(7,:)+polar_2nd(19,:)) & ! xxz xzx zxx
     &+DSIN(theta)**2*DCOS(theta)*DSIN(psi)*DCOS(psi) &
     & *(polar_2nd(6,:)+polar_2nd(8,:)+polar_2nd(12,:) & ! xyz xzy yxz
     & +polar_2nd(16,:)+polar_2nd(20,:)+polar_2nd(22,:)) & ! yzx zxy zyx
     &+DSIN(theta)*DSIN(psi) &
     & *(polar_2nd(14,:)+polar_2nd(2,:)-polar_2nd(18,:) & ! yyy xxy yzz
     & -polar_2nd(24,:)) & ! zyz
     &+DSIN(theta)*DCOS(psi) &
     & *(-polar_2nd(1,:)-polar_2nd(13,:)+polar_2nd(9,:) & ! xxx yyx xzz
     & +polar_2nd(21,:)) & ! zxz
     &+DSIN(theta)**3*DSIN(psi) &
     & *(-polar_2nd(2,:)-polar_2nd(4,:)-polar_2nd(10,:) & ! xxy xyx yxx
     & +polar_2nd(18,:)+polar_2nd(24,:)+polar_2nd(26,:)) & ! yzz zyz zzy
     &+DSIN(theta)**3*DCOS(psi) &
     & *(polar_2nd(5,:)+polar_2nd(11,:)+polar_2nd(13,:) & ! xyy yxy yyx
     & -polar_2nd(9,:)-polar_2nd(21,:)-polar_2nd(25,:)) & ! xzz zxz zzx
     &+DSIN(theta)**3*DSIN(psi)**3 &
     & *(-polar_2nd(14,:)+polar_2nd(2,:)+polar_2nd(4,:) & ! yyy xxy xyx
     & +polar_2nd(10,:)) & ! yxx
     &+DSIN(theta)**3*DCOS(psi)**3 &
     & *(polar_2nd(1,:)-polar_2nd(5,:)-polar_2nd(11,:) & ! xxx xyy yxy
     & -polar_2nd(13,:))  ! yyx
! chi_XYZ
      polar_2nd_ave(5,:)=num_den*0.5 &
     &*DCOS(theta)**2*(polar_2nd(6,:)-polar_2nd(12,:)) & ! xyz yxz
     &+DSIN(theta)**2*DSIN(psi)**2 &
     & *(polar_2nd(20,:)-polar_2nd(8,:)) & ! zxy xzy
     &+DSIN(theta)**2*DCOS(psi)**2 &
     & *(polar_2nd(16,:)-polar_2nd(22,:)) & ! yzx zyx
     &+DSIN(theta)**2*DSIN(psi)*DCOS(psi) &
     & *(polar_2nd(7,:)-polar_2nd(17,:)-polar_2nd(19,:) & ! xzx yzy zxx
     & +polar_2nd(23,:)) & ! zyy
     &+DSIN(theta)*DCOS(theta)*DSIN(psi) &
     & *(polar_2nd(5,:)-polar_2nd(9,:)-polar_2nd(11,:) & ! xyy xzz yxy
     & +polar_2nd(21,:)) & ! zxz
     &+DSIN(theta)*DCOS(theta)*DCOS(psi) &
     & *(-polar_2nd(4,:)+polar_2nd(10,:)-polar_2nd(18,:) & ! xyx yxx yzz
     & +polar_2nd(24,:)) ! zyz
! chi_XZY
      polar_2nd_ave(6,:)=num_den*0.5 &
     &*DCOS(theta)**2*(polar_2nd(8,:)-polar_2nd(16,:)) & ! xzy yzx
     &+DSIN(theta)**2*DSIN(psi)**2 &
     & *(polar_2nd(22,:)-polar_2nd(6,:)) & ! zyx xyz
     &+DSIN(theta)**2*DCOS(psi)**2 &
     & *(polar_2nd(12,:)-polar_2nd(20,:)) & ! yxz zxy
     &+DSIN(theta)**2*DSIN(psi)*DCOS(psi) &
     & *(polar_2nd(3,:)-polar_2nd(15,:)-polar_2nd(19,:) & ! xxz yyz zxx
     & +polar_2nd(23,:)) & ! zyy
     &+DSIN(theta)*DCOS(theta)*DSIN(psi) &
     & *(polar_2nd(5,:)-polar_2nd(9,:)-polar_2nd(13,:) & ! xyy xzz yyx
     & +polar_2nd(25,:)) & ! zzx
     &+DSIN(theta)*DCOS(theta)*DCOS(psi) &
     & *(-polar_2nd(2,:)+polar_2nd(10,:)-polar_2nd(18,:) & ! xxy yxx yzz
     & +polar_2nd(26,:))  ! zzy
!chi_ZXY
      polar_2nd_ave(7,:)=num_den*0.5 &
     &*DCOS(theta)**2*(polar_2nd(20,:)-polar_2nd(22,:)) & ! zxy zyx
     &+DSIN(theta)**2*DSIN(psi)**2 &
     & *(polar_2nd(16,:)-polar_2nd(12,:)) & ! yzx yxz
     &+DSIN(theta)**2*DCOS(psi)**2 &
     & *(polar_2nd(6,:)-polar_2nd(8,:)) & ! xyz xzy
     &+DSIN(theta)**2*DSIN(psi)*DCOS(psi) &
     & *(polar_2nd(3,:)-polar_2nd(7,:)-polar_2nd(15,:) & ! xxz xzx yyz
     & +polar_2nd(17,:)) & ! yzy
     &+DSIN(theta)*DCOS(theta)*DSIN(psi) &
     & *(polar_2nd(11,:)-polar_2nd(13,:)-polar_2nd(21,:) & ! yxy yyx zxz
     & +polar_2nd(25,:)) & ! zzx
     &+DSIN(theta)*DCOS(theta)*DCOS(psi) &
     & *(-polar_2nd(2,:)+polar_2nd(4,:)-polar_2nd(24,:) & ! xxy xyx zyz
     & +polar_2nd(26,:)) ! zzy
!
      DEALLOCATE(polar_2nd)
!
! Now calculate the experimentally measurable terms:
! chi_ssp, chi_sps chi_pss, chi_ppp, chi_spp, chi_pps, chi_psp
! Assuming Fresnel factors are independent of frequency of field
      ALLOCATE(chi_ssp(1:nstep))
      ALLOCATE(chi_sps(1:nstep))
      ALLOCATE(chi_pss(1:nstep))
      ALLOCATE(chi_ppp(1:nstep))
      ALLOCATE(chi_spp(1:nstep))
      ALLOCATE(chi_pps(1:nstep))
      ALLOCATE(chi_psp(1:nstep))
      fresnel=(/1.d0,1.d0,1.d0/)
      WRITE(*,*) "Frequency of the visible beam (cm-1):"
      READ(*,*) photon1
      WRITE(*,*) "Incident angles of the visible &
     & and IR beams (degree):"
      READ(*,*) beta1,beta2
! calculate beta by omega*sin(beta)=omega1*sin(beta1)+omega2*sin(beta2)
      beta1=beta1/PI*180.d0
      beta2=beta2/PI*180.d0
      DO i=1,nstep
        photon2=ibegin+(i-1)*step
        beta=DASIN((photon1*DSIN(beta1)+photon2*DSIN(beta2))/ &
     &(photon1+photon2))
!        WRITE(*,*) beta/PI*180.d0  
        chi_ssp(i)=fresnel(2)*fresnel(2)*fresnel(3) & !yyz
     &*DSIN(beta2)*polar_2nd_ave(4,i) !chi_YYZ
        chi_sps(i)=fresnel(2)*fresnel(3)*fresnel(2) & !yzy
     &*DSIN(beta1)*polar_2nd_ave(3,i) !chi_YZY
        chi_pss(i)=fresnel(3)*fresnel(2)*fresnel(2) & !zyy
     &*DSIN(beta)*polar_2nd_ave(2,i) !chi_ZYY
        chi_ppp(i)=-fresnel(1)*fresnel(1)*fresnel(3) & !xxz
     &*DCOS(beta)*DCOS(beta1)*DSIN(beta2)*polar_2nd_ave(4,i)  & !chi_XXZ
     &-fresnel(1)*fresnel(3)*fresnel(1) & !xzx
     &*DCOS(beta)*DSIN(beta1)*DCOS(beta2)*polar_2nd_ave(3,i)  & !chi_XZX
     &+fresnel(3)*fresnel(1)*fresnel(1) & !zxx
     &*DSIN(beta)*DCOS(beta1)*DCOS(beta2)*polar_2nd_ave(2,i)  & !chi_ZXX
     &+fresnel(3)*fresnel(3)*fresnel(3) & !zzz
     &*DSIN(beta)*DSIN(beta1)*DSIN(beta2)*polar_2nd_ave(1,i) !chi_ZZZ 
! chiral polarizaiton combination terms
        chi_spp(i)=fresnel(2)*fresnel(3)*fresnel(1) & !yzx
     &*DSIN(beta1)*DCOS(beta2)*(-polar_2nd_ave(6,i)) & !chi_YZX
     &+fresnel(2)*fresnel(1)*fresnel(3) & !yxz
     &*DCOS(beta1)*DSIN(beta2)*(-polar_2nd_ave(5,i))  !chi_YXZ
        chi_pps(i)=fresnel(3)*fresnel(1)*fresnel(2) & !zxy
     &*DSIN(beta)*DCOS(beta1)*polar_2nd_ave(7,i) & !chi_ZXY
     &-fresnel(1)*fresnel(3)*fresnel(2) & !xzy
     &*DCOS(beta)*DSIN(beta1)*polar_2nd_ave(6,i) !chi_XZY
        chi_psp(i)=fresnel(3)*fresnel(2)*fresnel(1) & !zyx
     &*DSIN(beta)*DCOS(beta2)*(-polar_2nd_ave(7,i)) & !chi_ZYX
     &-fresnel(1)*fresnel(2)*fresnel(3) & !xyz
     &*DCOS(beta)*DSIN(beta2)*polar_2nd_ave(5,i) !chi_XYZ
      ENDDO
!
      OPEN(UNIT=7,FILE="chi_exp",STATUS="UNKNOWN")
      WRITE(7,*) "Freqency/cm-1       chi_ssp              chi_ppp"
      DO i=1,nstep
        photon2=ibegin+(i-1)*step
        WRITE(7,"(F10.2,2E20.10)") photon2,ABS(chi_ssp(i)), &
     &ABS(chi_ppp(i))
      ENDDO
      CLOSE(7)
!
      DEALLOCATE(polar_2nd_ave)
      DEALLOCATE(chi_ssp)
      DEALLOCATE(chi_sps)
      DEALLOCATE(chi_pss)
      DEALLOCATE(chi_ppp)
      DEALLOCATE(chi_spp)
      DEALLOCATE(chi_pps)
      DEALLOCATE(chi_psp)
      END
