#!/bin/sh

# collecting data from Gaussian 09 D.01,
# checked by Zhitao Shen, 2021/10/18 

# changed to automatically generate input file of 'polarizability_2nd.f90' code
# The relevant vibrational modes are fixed in this script and should be checked carefully.
# Zhitao Shen, 2022/5/2

# read file name
echo -e "Enter file name: \c "
read FILE
#echo "$FILE"

#
#echo "Frequency (cm-1):"
grep "Frequencies ---" "$FILE" | awk \
     '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $3,$4,$5,$6,$7}' > freq.dat
#echo "Reduced mass (AMU):"
grep "Reduced masses ---" "$FILE" | awk \
     '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $4,$5,$6,$7,$8}' > reducedmass.dat
#echo "IR intensity (KM/Mole):"
#grep "IR Intensities ---" "$FILE" | awk \
#     '{printf "%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", $4,$5,$6,$7,$8}' > intensity.dat
#echo "Dipole derivative: "
grep "Dipole derivatives wrt mode" "$FILE" | awk '{printf "%15s%15s%15s\n", $6,$7,$8}'\
     > ddipole.dat
#echo "Polarizability derivative: "
grep -A4 "Polarizability derivatives wrt mode" "$FILE" \
     | grep -v "1             2             3" \
     | grep -v "\-\-" \
     | grep -v "Polarizability" \
     | awk '{printf "%15s%15s%15s\n", $2,$3,$4}' > dpolar.dat
#
sed -n -e '23, 25p' -e 30p freq.dat > temp1
sed -n -e '23, 25p' -e 30p reducedmass.dat > temp2
paste temp1 temp2 > temp
#
# prepare the input file for 'polarizability_2nd.f90' code
#
echo "Number of mode, begin and end of frequency, step size, damp(cm-1), ipolar(0 or 1), and ifresnel(0 or 1):" > temp.in
echo "4  2500  4000  1.0  3.0  0  0" >> temp.in
echo "Frequencies (cm-1) and reduced masses (amu):" >> temp.in
cat temp >> temp.in

rm freq.dat reducedmass.dat
rm temp1 temp2 temp
#
sed -n -e '23, 25p' -e 30p ddipole.dat > temp3
sed -n -e '67, 75p' -e '88, 90p' dpolar.dat > temp4 
#
echo "Dipole moment derivatives wrt normal modes (units is different in Gaussian [(km/mol)^1/2] and CP2K (au) output):" >> temp.in
echo "x          y        z" >> temp.in
cat temp3 >> temp.in
#
echo "Polarizability derivatives (A^2*amu^(-1/2)), order xx xy xz yx yy yz zx zy zz:" >> temp.in
echo " "
echo "******NOTICE******"
echo " "
echo "Polarizability derivatives (A^2*amu^(-1/2)), order xx xy xz yx yy yz zx zy zz:"
echo " "
echo "**********************"
echo " "
echo "Modes 23, 24, 25 and 30 are readin. Please check it!" >> temp.in
cat temp4 >> temp.in
#
rm ddipole.dat dpolar.dat
rm temp3 temp4
mv temp.in polarizability_2nd_g09.in 
#
echo 'Done! ^_^'
