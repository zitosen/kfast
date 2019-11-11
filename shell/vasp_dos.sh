#!/bin/bash

#########################################################################################
# Author  : zito
# Date    : 2018/10/24
# About   : Script to deal with the DOSCAR for VASP code. It's functions is similar with
#           the script 'split_dos', but they are different! Interactive is the most peculiar 
#           feature of the present script 'vasp_dos.sh', which ensure the data obtained are
#           just what you want.
# Version:  0.0
##########################################################################################

workdir=`pwd`
echo "Now you are in $workdir"

ls | grep DOSCAR > /dev/null
if [ $? == 0 ]
then
   echo "DOSCAR is found!"
else
   echo "There's no DOSCAR. Check it!"
fi
