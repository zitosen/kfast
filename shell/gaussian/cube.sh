#! /bin/bash

for filename in BDFF-T BDFT-T PBDFF-BT PBDFT-BT
do
   echo $filename'.fchk'
  cubegen 0 mo=homo $filename'.fchk' $filename'_homo.cube' -3  h
  cubegen 0 mo=lumo $filename'.fchk' $filename'_lumo.cube' -3  h
  cubegen 0 density=scf $filename'.fchk' $filename'_den.cube' -3  h
  cubegen 0 potential=scf $filename'.fchk' $filename'_pot.cube' -3  h
done 
