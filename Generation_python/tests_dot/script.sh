#!/bin/sh

for file in *.dot  
do
	if [ ! -f $file.pdf ]
	then
		dot -Tpdf $file -o $file.pdf
	fi
done
