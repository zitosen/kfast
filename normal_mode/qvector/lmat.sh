#!/bin/sh

echo -e "Enter file name: \c "
read FILE
echo -e "Enter dimension of l matrix: \c "
read L_DIM
BLOCK=`expr $L_DIM / 10 + 1`
#echo $BLOCK

LINE_1=`grep -n "EIGENVECTORS OF THE MASS-WEIGHTED" "$FILE" | cut  -d  ":"  -f  1`
LINE_1=`expr $LINE_1 + 1`
#echo $LINE_1
LINE_TOT=`expr 6 \* $BLOCK  + $L_DIM \* $BLOCK + \( $BLOCK - 1 \) \* $BLOCK + $LINE_1`
!echo $LINE_TOT

sed -n "$LINE_1,$LINE_TOT p" "$FILE" >> temp

sed '/^$/d' temp >> lmat
rm temp

