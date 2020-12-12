!#/bin/bash

head -n 1 ERR315857/species/species_profile.txt > top_species.txt
for file in $(ls)
do
	sed '1d' $file/species/species_profile.txt | head -n 1 >> top_species.txt
done
