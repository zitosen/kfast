      PROGRAM Euler_rotation
      IMPLICIT NONE
      INTEGER, PARAMETER       :: sr=4,sc=8,dr=8,dc=16
      REAL(KIND=dr), PARAMETER :: amu2au=1.82288848426455D3,    &
                                  au2cm=2.1947463137054D5,      &
                                  bohr2A=5.2917720859D-1, &
                                  PI=3.14159265358979323846264338328d0
      INTEGER                  :: natom,i,l,m,n,l_phi,m_theta,n_psi
      CHARACTER(LEN=100)       :: buffer
      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:,:) &
     &                         :: atom_label
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:,:) &
     &                         :: coor_xyz
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) &
     &                         :: coor_xyz_optimal
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:) &
     &                         :: atom_mass
      REAL(KIND=dr)            :: element(118,3),sum_coor_xyz(3,2), &
     & sum_atom_mass,center_of_mass(3,2),rot_euler(3,3),phi,theta,psi, &
     & phi_ini,phi_final,phi_step,theta_ini,theta_final,theta_step, &
     & psi_ini,psi_final,psi_step,temp,tot_square,phi_target, &
     & theta_target,psi_target
!    
! Assignment of atom mass in AMU, refer to 
! https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses                       !
      element=0.d0
      element(1,1)=1.0078250322d0
      element(6,1)=12.0000000000d0
      element(7,1)=14.0030740044d0
      element(8,1)=15.9949146196d0
      element(15,1)=30.9737619984d0
      element(16,1)=31.972071174d0
      element(34,1)=79.916522d0
! Read molecular orientation "geom1" on the surface, i.e., in XYZ frame
! Read molecular orientation "geom2" at xyz frame
      OPEN(UNIT=7,FILE="geom1",STATUS="OLD")
      OPEN(UNIT=77,FILE="geom2",STATUS="OLD")      
      READ(7,*) natom
      WRITE(*,*) "Atom number of geom1: ",natom
      READ(77,*) natom
      WRITE(*,*) "Atom number of geom2: ",natom
      ALLOCATE(atom_label(1:natom,1:2))
      ALLOCATE(atom_mass(1:natom))
      ALLOCATE(coor_xyz(1:3,1:natom,1:2))
      sum_coor_xyz=0.d0
      center_of_mass=0.d0
      sum_atom_mass=0.d0
      DO i=1,natom
        READ(7,*) atom_label(i,1),coor_xyz(1,i,1),coor_xyz(2,i,1),&
     &  coor_xyz(3,i,1)
        READ(77,*) atom_label(i,2),coor_xyz(1,i,2),coor_xyz(2,i,2),&
     &  coor_xyz(3,i,2)
        IF("H"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(1,1)
        ELSE IF("C"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(6,1)
        ELSE IF("N"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(7,1)
        ELSE IF("O"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(8,1)
        ELSE IF("P"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(15,1)
        ELSE IF("S"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(16,1)
        ELSE IF("Se"==TRIM(atom_label(i,1))) THEN
          atom_mass(i)=element(34,1)
        ELSE
          WRITE(*,*) "The atom mass of element ",atom_label(i,1), &
     &"is not assignment!"
        ENDIF
! calculate center of mass
        sum_coor_xyz(:,1)=sum_coor_xyz(:,1)+coor_xyz(:,i,1)*atom_mass(i)
        sum_coor_xyz(:,2)=sum_coor_xyz(:,2)+coor_xyz(:,i,2)*atom_mass(i)
        sum_atom_mass=sum_atom_mass+atom_mass(i)
      ENDDO
      center_of_mass(:,1)=sum_coor_xyz(:,1)/sum_atom_mass
      center_of_mass(:,2)=sum_coor_xyz(:,2)/sum_atom_mass
      WRITE(*,*) "Center of mass of geom1:"
      WRITE(*,*) center_of_mass(:,1)
      WRITE(*,*) "Center of mass of geom2:"
      WRITE(*,*) center_of_mass(:,2)
! check atom order
      DO i=1,natom
        IF(atom_label(i,1).NE.atom_label(i,2)) THEN
          WRITE(*,*) "ERROE! The atom order is different for geom1 and g&
     &eom2!"
          STOP
        ENDIF
      ENDDO
      CLOSE(7)
      CLOSE(77)
! New Cartesian coordinates setting center of mass as origin
      sum_coor_xyz=0.d0
      OPEN(UNIT=7,FILE="temp_geom.com",STATUS="UNKNOWN")
      WRITE(7,*) "#"
      WRITE(7,*) " "
      WRITE(7,*) "geom1 fixed at surface and rotated geom2"
      WRITE(7,*) " "
      WRITE(7,*) "1 1"
! geom1
      DO i=1,natom
        coor_xyz(:,i,1)=coor_xyz(:,i,1)-center_of_mass(:,1)
        sum_coor_xyz(:,1)=sum_coor_xyz(:,1)+coor_xyz(:,i,1)*atom_mass(i)
        WRITE(7,"(A4,3F18.10)") atom_label(i,1),coor_xyz(1,i,1), &
     &coor_xyz(2,i,1),coor_xyz(3,i,1)
      ENDDO
      center_of_mass(:,1)=sum_coor_xyz(:,1)/sum_atom_mass
      WRITE(*,*) "New center of mass of geom1:"
      WRITE(*,*) center_of_mass(:,1)
! geom2
      DO i=1,natom
        coor_xyz(:,i,2)=coor_xyz(:,i,2)-center_of_mass(:,2)
        sum_coor_xyz(:,2)=sum_coor_xyz(:,2)+coor_xyz(:,i,2)*atom_mass(i)
        WRITE(7,"(A4,3F18.10)") atom_label(i,1),coor_xyz(1,i,2), &
     &coor_xyz(2,i,2),coor_xyz(3,i,2)
      ENDDO
      center_of_mass(:,2)=sum_coor_xyz(:,2)/sum_atom_mass
      WRITE(*,*) "New center of mass of geom2:"
      WRITE(*,*) center_of_mass(:,2)
      WRITE(7,*) " "
      CLOSE(7)
!
! Maximum overlap of geom1(X,Y,Z) and geom2(x,y,z) by rotating geom2 
! through Euler matrix, 'rot_euler', that is, (X Y Z)^T=rot_euler(:,:)*(x y z)^T
! Euler angles, phi, theta and psi are defined with Z-y'-z convention
! refer to A. J. Moad and G. J. Simpson, J. Phys. Chem. B 2004, 108, 3548
      WRITE(*,*) "Initial value, final values and stepsize of &
     &phi (0-360 deg)"
      READ(*,*) phi_ini,phi_final,phi_step
      WRITE(*,*) "Initial value, final values and stepsize of &
     &theta (0-180 deg)"
      READ(*,*) theta_ini,theta_final,theta_step
      WRITE(*,*) "Initial value, final values and stepsize of &
     &psi (0-360 deg)"
      READ(*,*) psi_ini,psi_final,psi_step
      l_phi=INT((phi_final-phi_ini)/phi_step)
      m_theta=INT((theta_final-theta_ini)/theta_step)
      n_psi=INT((psi_final-psi_ini)/psi_step)
      ALLOCATE(coor_xyz_optimal(1:3,1:natom))
      coor_xyz_optimal=0.0d0
      tot_square=1.0d10
!      DO l=1,l_phi
!        phi=(phi_ini+(l-1)*phi_step)/180.d0*PI
!        DO m=1,m_theta
!          theta=(theta_ini+(m-1)*theta_step)/180.d0*PI
!          DO n=1,n_psi
!            psi=(psi_ini+(n-1)*psi_step)/180.d0*PI
            phi=phi_ini/180.d0*PI
            theta=theta_ini/180.d0*PI
            psi=psi_ini/180.d0*PI
            rot_euler(1,1)=-DSIN(psi)*DSIN(phi)+DCOS(theta)*DCOS(psi)*DCOS(phi)
            rot_euler(2,1)=DSIN(psi)*DCOS(phi)+DCOS(theta)*DCOS(psi)*DSIN(phi)
            rot_euler(3,1)=-DSIN(theta)*DCOS(psi)
            rot_euler(1,2)=-DCOS(psi)*DSIN(phi)-DCOS(theta)*DSIN(psi)*DCOS(phi)
            rot_euler(2,2)=DCOS(psi)*DCOS(phi)-DCOS(theta)*DSIN(psi)*DSIN(phi)
            rot_euler(3,2)=DSIN(theta)*DSIN(psi)
            rot_euler(1,3)=DSIN(theta)*DCOS(phi)
            rot_euler(2,3)=DSIN(theta)*DSIN(phi)
            rot_euler(3,3)=DCOS(theta)
            temp=0.d0
            DO i=1,natom
              coor_xyz(:,i,2)=MATMUL(rot_euler,coor_xyz(:,i,2))
              temp=temp+(coor_xyz(1,i,2)-coor_xyz(1,i,1))**2 &
     &                 +(coor_xyz(2,i,2)-coor_xyz(2,i,1))**2 &
     &                 +(coor_xyz(3,i,2)-coor_xyz(3,i,1))**2
            ENDDO
!            IF(temp.LT.tot_square) THEN
              tot_square=temp
              psi_target=psi/PI*180.d0
              theta_target=theta/PI*180.d0
              phi_target=phi/PI*180.d0
              coor_xyz_optimal(:,:)=coor_xyz(:,:,2)
!            ENDIF
!          ENDDO
!        ENDDO
!      ENDDO
      WRITE(*,*) 'The optimal phi, theta and psi are ', &
     &phi_target,theta_target,psi_target
      WRITE(*,*) 'Sum of distance square is ',tot_square
      OPEN(UNIT=7,FILE='final_geom.com',STATUS='UNKNOWN')
      WRITE(7,*) "#"
      WRITE(7,*) " "
      WRITE(7,*) "geom1 fixed at surface and rotated geom2"
      WRITE(7,*) " "
      WRITE(7,*) "1 1"
      DO i=1,natom
        WRITE(7,"(A4,3F18.10)") atom_label(i,1),coor_xyz(1,i,1), &
     &coor_xyz(2,i,1),coor_xyz(3,i,1)
      ENDDO
      DO i=1,natom
        WRITE(7,"(A4,3F18.10)") atom_label(i,1),coor_xyz_optimal(1,i), &
     &coor_xyz_optimal(2,i),coor_xyz_optimal(3,i)
      ENDDO
        WRITE(7,*) " "
      CLOSE(7)
      DEALLOCATE(atom_label)
      DEALLOCATE(atom_mass)
      DEALLOCATE(coor_xyz)
      DEALLOCATE(coor_xyz_optimal)
!
      END
