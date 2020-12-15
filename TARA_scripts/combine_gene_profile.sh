!#/bin/bash

head -n 1 ERR315857/genes/summary.txt > genes_merge.txt
for file in $(ls)
do
	sed '1d' $file/genes/summary.txt | head -n 50 >> genes_merge.txt 
done
