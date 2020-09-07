#!/bin/sh
####################################################################
#
# define the following variables according to your needs
#
  outdir=./tmp
  pseudo_dir=../pseudo
# espresso_dir=/where your espresso is installed/
####################################################################

if [ -f si.etot_vs_ecut ] ; then
   rm si.etot_vs_ecut
fi
if [ -f si.ecutwfc.in ] ; then
   rm si.ecutwfc.in
fi
if [ -f si.ecutwfc.out ] ; then
   rm si.ecutwfc.out
fi
if [ -d tmp ] ; then
   rm -rf tmp
fi

for ecutwfc in 12.0 16.0 20.0 24.0 28.0 32.0 36.0 40.0; do

# self-consistent calculation
cat > si.ecutwfc.in << EOF
 &control
    prefix='silicon',
    pseudo_dir = '$pseudo_dir/',
    outdir='$outdir/'
 /
 &system    
    ibrav=  2, celldm(1) =10.2, nat=  2, ntyp= 1,
    ecutwfc = $ecutwfc, 
 /
 &electrons
 /
ATOMIC_SPECIES
 Si  28.086  Si.pz-vbc.UPF
ATOMIC_POSITIONS
 Si 0.00 0.00 0.00 
 Si 0.25 0.25 0.25 
K_POINTS automatic
   2 2 2 1 1 1
EOF

# If pw.x is not found, specify the correct value for $espresso_dir,
# use $espresso_dir/bin/pw.x instead of pw.x

mpirun pw.x -in si.ecutwfc.in > si.ecutwfc.out

#grep -e 'lattice parameter' -e ! si.eos.out | \
#      awk '/lattice/{alat=$(NF-1)}/!/{print alat, $(NF-1)}' >> si.etot_vs_alat
grep -E '!.*total energy'  si.ecutwfc.out | \
       awk -v ecut="$ecutwfc" '{print ecut, "\t" , $(NF-1)}' >> si.etot_vs_ecut
done
