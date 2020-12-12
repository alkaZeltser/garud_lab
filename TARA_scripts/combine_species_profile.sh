!#/bin/bash

head -n 1 ERR315857/species/species_profile.txt > species_merge.txt
for file in $(ls)
do
	sed '1d' $file/species/species_profile.txt | head -n 300 >> species_merge.txt 
done
