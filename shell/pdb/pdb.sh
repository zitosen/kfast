#! /bin/sh

# prepare input file for QM/MM calculations
# by zito  2018/12/29

# read file names
echo -n "The first file: "
read filename1
echo -n "The second file: "
read filename2

# output file names
fileout=out.pdb

# handle main data
sed '/END/d' $filename1 | sed '/REMARK/d' > tmp1
sed '/END/d' $filename2 | sed '/REMARK/d' > tmp2
cat tmp1 tmp2 >> tmp3
nl tmp3 | awk '{$3 = $1; printf "%4s %6s %4s %3s %5s %11s %7s %7s %5s %5s %11s\n", \
                                 $2, $3, $4, $5, $7, $8, $9, $10, $11, $12, $13}' > tmp4
# add head and tail
head -n 2 $filename1 > tmp_head
tail -n 1 $filename1 > tmp_tail
cat tmp_head tmp4 tmp_tail >> tmp5
sed 's///' tmp5 > $fileout

# remove temporary files
rm tmp*

echo "Done!"
