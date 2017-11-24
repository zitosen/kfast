!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE ir_raman_tensor
      IMPLICIT NONE
      SAVE
      INTEGER :: i_resonant
      REAL(KIND=8) :: omega(3)
      CHARACTER(100) :: file_TDM,file_E
!      
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     IR and Raman active tensor calculations                           !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE tensor_cal(alpha,nstate)
      IMPLICIT NONE
      REAL(KIND=8),INTENT(OUT) :: alpha(1:3,1:3)
      INTEGER,INTENT(IN) :: nstate
      REAL(KIND=8) :: mu(0:nstate,0:nstate,1:3),energy(0:nstate)
      INTEGER :: istate,jstate,i,j
      LOGICAL :: ok
!
      CALL dipole_moment(mu,nstate)      !mu in au
      CALL energy_cis_d(energy,nstate)   !energy in au
!
!!!!! Off-resonance-resonance case !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      alpha=0.d0
      DO i=1,3                    !x,y,z
        DO j=1,3                  !x,y,z
          DO istate=0,nstate
            alpha(i,j)= alpha(i,j)+ &
     &                 (mu(i_resonant,istate,j)*mu(istate,0,i))/ &
     &                 (energy(istate)-omega(1))+ &
     &                 (mu(i_resonant,istate,i)*mu(istate,0,j))/ &
     &                 (energy(istate)-omega(2))
          END DO
        END DO
      END DO
      END SUBROUTINE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Subroutine 'dipole_moment' reads dipole moment from log file of       !
! Gaussian (09) package.                                                !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dipole_moment(mu,nstate)
      IMPLICIT NONE
      REAL(KIND=8),INTENT(OUT) :: mu(0:nstate,0:nstate,1:3)
      INTEGER,INTENT(IN) :: nstate
      REAL(KIND=8) :: temp1,temp2
      CHARACTER(100) :: mu_mg_flag,mu_mn_flag,buffer
      INTEGER :: istate,jstate
      LOGICAL :: ok
!
      mu=0.d0                                                        !au
      mu(0,0,1:3)=(/0.d0,0.d0,-1.227105d0/)                          !au
      mu(i_resonant,i_resonant,1:3)=(/0.d0,0.d0,-0.669263d0/)        !au
      mu(0,i_resonant,1:3)=(/0.d0,0.d0,0.d0/)                        !au
      mu(i_resonant,0,1:3)=(/0.d0,0.d0,0.d0/)                        !au
!
      mu_mg_flag=" Ground to excited state transition electric dipole mo&
     &ments (Au):"
      mu_mn_flag=" Excited to excited state transition electric dipole m&
     &oments (Au):"
      INQUIRE(file=file_TDM,exist=ok)
      IF(ok) THEN
        OPEN(unit=77,file=file_TDM,status='old',form='formatted',       &
     &  access='sequential',blank='null')
        DO WHILE(.TRUE.)
          READ(77,'(A100)') buffer
          IF(mu_mg_flag.EQ.TRIM(buffer)) THEN
            READ(77,'(A100)') buffer
            istate=1
            DO WHILE(istate.LT.nstate)
              READ(77,*) istate,mu(istate,0,1),mu(istate,0,2),          &
     &        mu(istate,0,3),temp1,temp2      !0 is ground state
              mu(0,istate,1)=mu(istate,0,1)
              mu(0,istate,2)=mu(istate,0,2)
              mu(0,istate,3)=mu(istate,0,3)
            END DO
          END IF
          IF(mu_mn_flag.EQ.TRIM(buffer)) THEN
            READ(77,'(A100)') buffer
            istate=2
            jstate=1
            DO WHILE(jstate.LT.(nstate-1))
              READ(77,*) istate,jstate,mu(istate,jstate,1),             &
     &        mu(istate,jstate,2),mu(istate,jstate,3),temp1,temp2
              mu(jstate,istate,1)=mu(istate,jstate,1)
              mu(jstate,istate,2)=mu(istate,jstate,2)
              mu(jstate,istate,3)=mu(istate,jstate,3)
            END DO
            EXIT
          END IF
        END DO
        CLOSE(77)
      ELSE
        WRITE(*,*) "File ",'"',TRIM(file_TDM),'"', " doesn't exist."
      END IF
      RETURN
      END SUBROUTINE dipole_moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Subroutine 'energy_cis_d' reads energies of CIS(D) calculations       !
! by Gaussian (09) package.                                             !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE energy_cis_d(energy,nstate)
      IMPLICIT NONE
      REAL(KIND=8),INTENT(OUT) :: energy(0:nstate)
      INTEGER,INTENT(IN) :: nstate
      INTEGER :: i_au,istate
      CHARACTER(LEN=100) :: cis_d_flag,buffer
      LOGICAL :: ok
!
      energy(0)=0.d0     !au
      INQUIRE(file=file_E,exist=ok)
      IF(ok) THEN
        cis_d_flag=" CIS(D) Exc. E:"
        OPEN(unit=78,file=file_E,status='old',form='formatted',access=  &
     &  'sequential',blank='null')
        istate=1
        DO WHILE(istate.LE.nstate)
          READ(78,'(A100)') buffer
          IF (cis_d_flag.EQ.buffer(1:15)) THEN
            i_au=INDEX(buffer,"a.u.")
            READ(buffer(i_au-21:i_au),'(F20.11)') energy(istate)
!            WRITE(*,'(F18.11,2X,A4)') energy(istate),"a.u."
            istate=istate+1
          END IF
        END DO
        CLOSE(78)
      ELSE
        WRITE(*,*) "File ",'"',TRIM(file_E),'"'," doesn't exist."
      END IF
      RETURN
      END SUBROUTINE energy_cis_d
      END MODULE ir_raman_tensor
