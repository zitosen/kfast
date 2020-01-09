#!/bin/bash
#
#########################################################################
#
#     Objetive: Extracting data from Gaussian 09 output
#               Note keywords "freq(raman,hpmodes) iop(7/33=1,7/33=3)" 
#               should be specifed in the route section in Gaussian 
#               calculation
#                               
#     Author  : Z. Shen
#     Email   : shenzt@vip.henu.edu.cn
#     Version : 0.0
#     Date    : 2020/1/2
#
#########################################################################
#
# read file name
  echo -e "Enter file name: \c "
  read FILE
#echo "$FILE"

#
  echo "Frequency (cm-1):"
  grep "Frequencies ---" "$FILE" | awk \
       '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $3,$4,$5,$6,$7}' > freq.dat
  echo "Reduced mass (AMU):"
  grep "Reduced masses ---" "$FILE" | awk \
       '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $4,$5,$6,$7,$8}' > reducedmass.dat
  echo "IR intensity (KM/Mole):"
  grep "IR Intensities ---" "$FILE" | awk \
       '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $4,$5,$6,$7,$8}' > intensity.dat
  echo "Dipole derivative: "
  grep "Dipole derivatives wrt mode" "$FILE" | awk '{printf "%15s%15s%15s\n", $6,$7,$8}'\
       > ddipole.dat
  echo "Polarizability derivative: "
  grep -A4 "Polarizability derivatives wrt mode" "$FILE" \
       | grep -v "1             2             3" \
       | grep -v "\-\-" \
       | grep -v "Polarizability" \
       | awk '{printf "%15s%15s%15s\n", $2,$3,$4}' > dpolar.dat
#
  echo "Extracting molecular geometry in Cartesian coordinates:"
  line_begin=`grep -n "Input orientation:" $FILE | sed 's/[^0-9]//g'`
  line_begin=`expr $line_begin + 5`
  line_end=`grep -n -A1000 "Input orientation:" $FILE | \
               grep '\-\-\-' | sed "s/[^0-9]//g" | sed -n '3,3p'`
  line_end=`expr $line_end - 1`
  sed -n "$line_begin,$line_end p" $FILE > geom_s0.dat 
#
  echo "DONE!"

 
