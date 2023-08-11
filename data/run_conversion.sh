#!/bin/bash 
set -e

# for n in 10 100 1000 10000 100000 1000000; do
for n in 1000; do
	echo "Subsetting for $n"

	stem="chr21_n$n"
	python3 ../src/convert.py subset-trees simulated_chrom_21.ts $n $stem.ts
	make $stem.sgz
done

