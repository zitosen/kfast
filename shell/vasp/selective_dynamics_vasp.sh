#!/bin/sh
# processing vasp POSCAR file for selective dynamics calculations
# by zito, shenzt@vip.henu.edu.cn, 2020/9/27

# atoms with z coordinate <  z_value will be fiexed
z_value=0.1

if [ -f temp ] ; then
   rm temp
fi

awk 'NR<10{print $0}' POSCAR > temp
awk -v z=$z_value 'NR>9{ if($3>z) printf "%15.9f %15.9f %15.9f  T  T  T\n",  $1, $2, $3 ; \
                   else printf "%15.9f %15.9f %15.9f  F  F  F\n",  $1, $2, $3 }' \
                   POSCAR >> temp
#
echo "Done!"
