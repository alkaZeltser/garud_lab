setwd("/Users/nikaz/Downloads")

library(tidyverse)
library(readxl)

#Read in patric metadata
patric = read_excel(path="PATRIC_metadata.xls", col_names=T)
patric = patric[69:714, ]

#read in Pachiadaki supplemental table S2
tropics = read_excel(path="gorg-tropics_sags_tableS2.xlsx", col_names=T)
tropics = subset(tropics, select = c(SAG, `Genome completeness (%)`))
colnames(tropics) = c("Strain", "Completeness")
colnames(tropics)

#left outer join to match up SAR86 patric genomes with their completeness.
master = merge(patric, tropics, all.x = T, sort = F)

#Filter for completeness over 50%
complete = subset(master, Completeness >= 50, select = c(`Genome ID`, Strain, Completeness))

#Create path file for fastANI
path_file = paste("/u/home/n/nzeltser/project-ngarud/sar86/", complete$`Genome ID`, ".fna", sep = "")

#Export path file
write.table(path_file, file = "complete_genomes.txt", row.names = F, col.names = F, quote = F)

## After fastANI analysis on hoffman ##
#Read in fastANI results
fastani = read.table(file = "fastani.out", as.is = T)
colnames(fastani) = c("Query", "Reference", "ANI", "mapped_fragment_count", "query_fragment_count")

#Remove redundant path section to make the data more human readable
fastani$Query = gsub("/u/home/n/nzeltser/project-ngarud/sar86/", "", fastani$Query)
fastani$Reference = gsub("/u/home/n/nzeltser/project-ngarud/sar86/", "", fastani$Reference)

## Filter out self pairings ##
no_self_pairs = subset(fastani, Query != Reference)

## Within Species Pairs (ANI > 95)
same_species = subset(no_self_pairs, ANI > 95)

## Trying to figure out how to weed out duplicate pairs.. not working
reorder = no_self_pairs[,c(2,1,3,4,5)]
colnames(reorder)[c(1,2)] = c("Query", "Reference")
test = rbind(no_self_pairs, reorder)
test = test[order(test$Query, test$Reference),]


## Graphing ANI distribution ##
# ANI distribution with only self pairings removed:

hist(no_self_pairs$ANI, plot = T, col = "tomato", breaks = 50, freq = T,
     main = "ANI Distribution", xlab = "Pairwise ANI `(%)`")


# \/u\/home\/n\/nzeltser\/project-ngarud\/sar86\/
