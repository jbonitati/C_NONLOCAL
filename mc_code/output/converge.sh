#!/bin/bash
EXPECT=0.60131267
outfile='convtest.txt'
rm -f $outfile
for N in {1..50}
do
  sed -i "/Basis_size/ c\
Basis_size = $N" config.ini
  OUTPUT=$(./run | tail -1)
  DIFF=$(echo "($EXPECT - $OUTPUT) / $EXPECT" | bc -l)
  echo "$N $DIFF">>$outfile
done
