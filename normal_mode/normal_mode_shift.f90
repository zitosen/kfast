!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate molecular displaced geometry along normal mode vibrational direction.
! The molecular equilibrium geometry is obtained by CP2K code.
!                     2016.10.5, by zito
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      PROGRAM NORMAL_MODE_SHIFT
      IMPLICIT NONE
      INTEGER, PARAMETER               :: dr=8,dc=16
      REAL(KIND=dr), PARAMETER         :: amu2au=1.82288848426455D3,    &
                                          au2cm=2.1947463137054D5,      &
                                          bohr2A=5.2917720859D-1
!      REAL(KIND=dr), PARAMETER         :: DQ=5.0D-3    ! unit is angstrom*amu^(1/2)
      REAL(KIND=dr), PARAMETER         :: DQ=1.0D-2
      REAL(KIND=dr), ALLOCATABLE, &
      DIMENSION(:)                     :: freq,red_mass
      REAL(KIND=dr), ALLOCATABLE, DIMENSION(:,:) & 
                                       :: coor_eq,coor_p,coor_m,lmatrix
      INTEGER                          :: nmode,natom,i,j
      CHARACTER(60)                    :: temp
      CHARACTER(20)                     :: mode_num,filename
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) &
                                       :: atom_label
!
      OPEN(7,FILE='normal_mode_shift.inp',STATUS='OLD')
      READ(7,*) nmode,natom
      ALLOCATE(freq(nmode))
      ALLOCATE(red_mass(nmode))
      ALLOCATE(atom_label(natom))
      ALLOCATE(coor_eq(natom,3))
      ALLOCATE(coor_p(natom,3))
      ALLOCATE(coor_m(natom,3))
      ALLOCATE(lmatrix(natom,3))
      READ(7,*) temp
! read in frequencies (cm-1) and reduced masses (amu)
      DO i=1,nmode
        READ(7,*) freq(i),red_mass(i)
      ENDDO
      READ(7,*) temp
! read in Cartestian coordinates at equilibrium geometry configuration
      DO i=1,natom
        READ(7,*) atom_label(i),coor_eq(i,1),coor_eq(i,2),coor_eq(i,3)
! unit transform, bohr-->angstrom
        coor_eq(i,:)=coor_eq(i,:)*bohr2A
      ENDDO
!
      READ(7,*) temp
! read in l matrix of normal mode
      DO j=1,nmode
        WRITE(mode_num,'(I4)') j
        filename='mode'//TRIM(ADJUSTL(mode_num))//'.p'
        OPEN(77,FILE=filename,STATUS='NEW')
        filename='mode'//TRIM(ADJUSTL(mode_num))//'.m'
        OPEN(78,FILE=filename,STATUS='NEW')
        lmatrix(:,:)=0.D0
        coor_p(:,:)=0.D0
        coor_m(:,:)=0.D0
        READ(7,*) temp
        DO i=1,natom
! here l matrix is dimensionless
          READ(7,*) lmatrix(i,1),lmatrix(i,2),lmatrix(i,3)
! calculate DX by the relation: DX = l * M^(-1/2) * DQ
! DX: displacment of Cartestian coordinate
! DQ: displacment of normal mode coordinate
! M: diagonal matrix, M_ii is the reduced mass of the i-th mode
          lmatrix(i,:)=lmatrix(i,:)/DSQRT(red_mass(j))
          coor_p(i,:)=coor_eq(i,:)+lmatrix(i,:)*DQ
          coor_m(i,:)=coor_eq(i,:)-lmatrix(i,:)*DQ
          WRITE(77,'(A10,3F18.8,3X)') atom_label(i),coor_p(i,:)
          WRITE(78,'(A10,3F18.8,3X)') atom_label(i),coor_m(i,:)
        ENDDO
        CLOSE(77)
        CLOSE(78)
      ENDDO
      CLOSE(7)
      DEALLOCATE(freq)
      DEALLOCATE(red_mass)
      DEALLOCATE(atom_label)
      DEALLOCATE(coor_eq)
      DEALLOCATE(coor_p)
      DEALLOCATE(coor_m)
      DEALLOCATE(lmatrix)
!
      END PROGRAM NORMAL_MODE_SHIFT
