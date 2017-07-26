#!/bin/bash
EXPECT=
outfile='convtest.txt'
rm -f $outfile
for N in {1..20}
do
  echo $N>>$outfile
  sed -i "/Basis_size/ c\
Basis_size = $N" config.ini
  OUTPUT=$(./run | tail -1)
  DIFF=($OUTPUT - $EXPECT)/($EXPECT)
  echo $DIFF>>$outfile
  echo "\n">>$outfile
done
