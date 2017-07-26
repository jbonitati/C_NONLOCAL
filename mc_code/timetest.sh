#!/bin/bash
outfile='timetest.txt'
rm -f $outfile
for N in 5 10 20 40 80 160 320 640 1280 2000 3000 4000 5000
do
  echo $N>>$outfile
  sed -i.bak "/Basis_size/ c\
Basis_size = $N" config.ini
  ./run | grep time >> "$outfile"
done
