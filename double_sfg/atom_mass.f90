!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!     atom_mass: assignment of atomic mass                              !
!                                                                       !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE atom_mass
      IMPLICIT NONE
!
      CONTAINS
      FUNCTION n2m(iso,num)
      IMPLICIT NONE
      REAL(KIND=8) :: amu(0:2,1:10),n2m
      INTEGER,INTENT(IN) :: iso,num !!iso: isotope type!
!
      n2m=0.d0
      amu(0,1:10)=(/1.007825d0,4.00260d0,7.016003d0,9.012182d0,&
     &11.009305d0,12.0d0,14.003074d0,15.994915d0,18.9984032d0,&
     &19.992435d0/)
      IF(0<iso<=2.AND.0<num<=10) THEN
        n2m=amu(iso,num)
      ELSE
        WRITE(*,'(A30)') 'Check the "atom_mass" module!'
      END IF
      END FUNCTION n2m
      END MODULE atom_mass

