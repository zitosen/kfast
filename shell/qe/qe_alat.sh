#!/bin/sh
####################################################################
#
# define the following variables according to your needs
#
  outdir=./tmp
  pseudo_dir=../pseudo
# espresso_dir=/where your espresso is installed/
####################################################################

if [ -f si.etot_vs_alat ] ; then
   rm si.etot_vs_alat
fi
if [ -f si.alat.in ] ; then
   rm si.alat.in
fi
if [ -f si.alat.out ] ; then
   rm si.alat.out
fi
if [ -d tmp ] ; then
   rm -rf tmp
fi

for alat in 9.8 9.9 10.0 10.1 10.2 10.3 10.4 10.5 10.6 10.7 10.8; do

# self-consistent calculation
cat > si.alat.in << EOF
 &control
    prefix='silicon',
    pseudo_dir = '$pseudo_dir/',
    outdir='$outdir/'
 /
 &system    
    ibrav=  2, celldm(1) =$alat, nat=  2, ntyp= 1,
    ecutwfc = 20, 
 /
 &electrons
 /
ATOMIC_SPECIES
 Si  28.086  Si.pz-vbc.UPF
ATOMIC_POSITIONS
 Si 0.00 0.00 0.00 
 Si 0.25 0.25 0.25 
K_POINTS automatic
   4 4 4  1 1 1
EOF

# If pw.x is not found, specify the correct value for $espresso_dir,
# use $espresso_dir/bin/pw.x instead of pw.x

mpirun -n 24 pw.x -in si.alat.in > si.alat.out

grep -E '!.*total energy'  si.alat.out | \
       awk -v alat="$alat" '{print alat, "\t" , $(NF-1)}' >> si.etot_vs_alat
done
