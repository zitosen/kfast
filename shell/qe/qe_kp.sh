#!/bin/sh
####################################################################
#
# define the following variables according to your needs
#
  outdir=./tmp
  pseudo_dir=../pseudo
# espresso_dir=/where your espresso is installed/
####################################################################

if [ -f si.etot_vs_nkp ] ; then
   rm si.etot_vs_nkp
fi
if [ -f si.kp.in ] ; then
   rm si.kp.in
fi
if [ -f si.kp.out ] ; then
   rm si.kp.out
fi
if [ -d tmp ] ; then
   rm -rf tmp
fi

for nkp in 2 3 4 5 6 7 8; do

# self-consistent calculation
cat > si.kp.in << EOF
 &control
    prefix='silicon',
    pseudo_dir = '$pseudo_dir/',
    outdir='$outdir/'
 /
 &system    
    ibrav=  2, celldm(1) =10.2, nat=  2, ntyp= 1,
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
   $nkp $nkp $nkp 1 1 1
EOF

# If pw.x is not found, specify the correct value for $espresso_dir,
# use $espresso_dir/bin/pw.x instead of pw.x

mpirun -n 24 pw.x -in si.kp.in > si.kp.out

grep -E '!.*total energy'  si.kp.out | \
       awk -v nkp="$nkp" '{print nkp, "\t" , $(NF-1)}' >> si.etot_vs_nkp
done
