setwd("/Users/nikaz/Desktop")

#Import SRA accessions from TARA project
project_sra_accessions = read.table(file = "TARA_metadata/TARA_accessions.txt", header = F, as.is = T)

#Import keys for SRA accession -> TARA alias conversion from ENA
key87 = read.table(file = "Tara_metadata/filereport_read_run_PRJEB1787_tsv.txt", sep = "\t", header = T, as.is = T)
key88 = read.table(file = "Tara_metadata/filereport_read_run_PRJEB1788_tsv.txt", sep = "\t", header = T, as.is = T)

#Join the two keys (different sample categories) into one key
key = rbind(key87, key88)

#Filter key on just the relevant sra accessions:
key = subset(key, run_accession %in% project_sra_accessions$V1)

#Import TARA metadata from Pangaea (note: had to manually edit file to delete data description)
metadata = read.table(file = "TARA_metadata/TARA_sample_enviro.tab", sep = "\t", header = F, as.is = T, fill = T)

#R was having trouble parsing the header, so I'm importing it seperately
meta_header = read.table(file = "TARA_metadata/TARA_header.txt", sep = "\t", header = F, as.is = T, 
                         comment.char = "", quote = "")
colnames(metadata) = meta_header[1,]

#Filter Metadata by relevant sample aliases (from key)
metadata = subset(metadata, `Sample ID (TARA_barcode#)` %in% key$sample_alias)

#Filter Metadata by variables related to geographic region 
# (descriptions can be found here: https://doi.pangaea.de/10.1594/PANGAEA.853810?format=html#download)

geo_data = metadata[,c(1,4,5,6,8,9,10,11,12,13,14,16,17,18,19)]

#Merge with sra accessions

geo_data = merge(key[,c(4,7)], geo_data, by.x = "sample_alias", by.y = "Sample ID (TARA_barcode#)")

#Write to file
write.table(geo_data, file = "TARA_metadata/custom_tara_geo_metadata.txt", quote = F, sep = "\t", row.names = F)
