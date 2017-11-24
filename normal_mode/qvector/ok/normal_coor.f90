      PROGRAM normal_coor
      IMPLICIT NONE
      CHARACTER(100) :: filename,buffer
      CHARACTER(40) :: fname,fnum,temp
      INTEGER :: i,j,m,n,imax
      INTEGER,PARAMETER :: lmat_dim=72,num_file=100
      REAL(KIND=8),ALLOCATABLE :: lmat(:,:),xyz(:,:),qnorm(:),mass(:),    &
                                  xyz_m(:,:)

      ALLOCATE(lmat(1:lmat_dim,1:lmat_dim))
      ALLOCATE(xyz(1:lmat_dim,1:num_file))
      ALLOCATE(qnorm(1:lmat_dim))
      ALLOCATE(mass(1:lmat_dim/3))
      ALLOCATE(xyz_m(1:lmat_dim,1:num_file))
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
      CLOSE(79)
!
      OPEN(unit=77,file="qvector",form='formatted',status='new')
      DO m=1,100
        WRITE(fnum,'(I4)') m
        fname="mndogeom"//TRIM(ADJUSTL(fnum))
        WRITE(*,*) fname
        OPEN(unit=79,file=fname,form='formatted',status='old')
          DO i=1,lmat_dim/3
            READ(79,*) mass(i),xyz((i-1)*3+1,m),temp,xyz((i-1)*3+2,m),temp, &
     &                 xyz((i-1)*3+3,m),temp
            xyz_m((i-1)*3+1,m)=xyz((i-1)*3+1,m)*mass(i)
            xyz_m((i-1)*3+2,m)=xyz((i-1)*3+2,m)*mass(i)
            xyz_m((i-1)*3+3,m)=xyz((i-1)*3+3,m)*mass(i)
          ENDDO
        CLOSE(79)
! calculate normal mode coordinates
        qnorm=MATMUL(TRANSPOSE(lmat),xyz_m(:,m))
        WRITE(77,'(A18,I4)') "For file number: ", m
        WRITE(77,'(F12.6)') qnorm
        WRITE(77,'(A)') " "
      ENDDO
      CLOSE(77)
!
      DEALLOCATE(lmat)
      DEALLOCATE(xyz)
      DEALLOCATE(qnorm)
      DEALLOCATE(mass)
      DEALLOCATE(xyz_m)
!
      END
