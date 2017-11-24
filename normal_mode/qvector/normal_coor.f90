! Extracting normal mode coordinates from MNDO2005 program!
! by zito
!
      PROGRAM normal_coor
      IMPLICIT NONE
      CHARACTER(100) :: filename,buffer,temp
      INTEGER :: i,j,m,n,imax
      INTEGER,PARAMETER :: lmat_dim=72
      REAL(KIND=8),ALLOCATABLE :: lmat(:,:),xyz(:),qnorm(:),mass(:),    &
                                  xyz_m(:)

      ALLOCATE(lmat(1:lmat_dim,1:lmat_dim))
      ALLOCATE(xyz(1:lmat_dim))
      ALLOCATE(qnorm(1:lmat_dim))
      ALLOCATE(mass(1:lmat_dim/3))
      ALLOCATE(xyz_m(1:lmat_dim))
      imax=lmat_dim/10+1
      filename="lmat"
      OPEN(unit=79,file=filename,form='formatted',status='old')
      DO i=1,imax
        IF(i.LT.imax) THEN
          READ(unit=79,fmt='(A100)') buffer
          READ(unit=79,fmt='(A100)') buffer
          DO j=1,lmat_dim
            READ(79,*)                                                  &
     &      temp,lmat(j,((i-1)*10+1)),lmat(j,((i-1)*10+2)),             &
     &lmat(j,((i-1)*10+3)),lmat(j,((i-1)*10+4)),lmat(j,((i-1)*10+5)),   &
     &lmat(j,((i-1)*10+6)),lmat(j,((i-1)*10+7)),lmat(j,((i-1)*10+8)),   &
     &lmat(j,((i-1)*10+9)),lmat(j,((i-1)*10+10))
          ENDDO
        ELSE
          READ(unit=79,fmt='(A100)') buffer
          READ(unit=79,fmt='(A100)') buffer
          DO j=1,lmat_dim
            READ(79,*) temp,lmat(j,((imax-1)*10+1)),lmat(j,((imax-1)*10+2))
          ENDDO
        ENDIF
      ENDDO
!      WRITE(*,*) lmat(72,71:72)
      CLOSE(79)
!
      OPEN(unit=79,file="xyz",form='formatted',status='old')
        DO i=1,lmat_dim/3
          READ(79,*) temp,mass(i),xyz((i-1)*3+1),xyz((i-1)*3+2),xyz((i-1)*3+3)
          xyz_m((i-1)*3+1)=xyz((i-1)*3+1)*mass(i)
          xyz_m((i-1)*3+2)=xyz((i-1)*3+2)*mass(i)
          xyz_m((i-1)*3+3)=xyz((i-1)*3+3)*mass(i)
        ENDDO
      CLOSE(79)
!
      qnorm=MATMUL(TRANSPOSE(lmat),xyz_m)
      OPEN(unit=79,file="qvector",form='formatted',status='new')
        WRITE(79,'(F12.6)') qnorm
      CLOSE(79)
!
      DEALLOCATE(lmat)
      DEALLOCATE(xyz)
      DEALLOCATE(qnorm)
      DEALLOCATE(mass)
      DEALLOCATE(xyz_m)
!
      END
