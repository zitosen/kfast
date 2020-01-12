!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!  References: (1)P. J. Mohr, B. N. Taylor, and D. B. Newell, “CODATA
!  Recommended Values of the Fundamental Physical Constants: 2010,” Rev.
!  Mod. Phys., 84 (2012) 1527-1605.
!                (2) P. J. Mohr, B. N. Taylor, and D. B. Newell, “CODATA
!  Recommended Values of the Fundamental Physical Constants: 2010,” Chem.
!  Ref. Data, 41 (2012) 043109.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE constants
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: AU2eV=27.2114d0,AMU2kg=1.660538921d-27,&
     &elec2kg=0.910938291d-30,AU2cm=219474.63d0,AU2kcal=627.5095d0,     &
     &Bohr2A=0.52917721092d0,AMU2elec=1822.8885d0,                      &
     &k_Boltzman=1.3806488d-23
      REAL(KIND=8), PARAMETER :: h_Planck=6.62606957d-34,               &
     &ONEMOLE=6.02214129d23,PI=3.14159265358979323846264338328d0,       &
     &clight=2.99792458d8
      END MODULE
