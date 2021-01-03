
setwd("/u/scratch/n/nzeltser/midas/snps")

#install.packages("ggfortify")

library(ggplot2)
library(ggfortify)

# Read in file containing name of every species 
species = read.table(file = "species.txt", header = T, as.is = T)

# Open graphing device
pdf("PCA_bgregions.pdf")

# Iterate over each species
for (s in species$species_name){
  
  #Import snp frequency file
  snp_freqs = read.table(file = paste(s, "/snps_ref_freq.txt", sep = ""), header = T, as.is = T)
  
  #Transpose data
  flipped_snps = data.frame(t(snp_freqs[-1]))
  colnames(flipped_snps) = snp_freqs[,1]
  
  #Import metadata
  geo_data = read.table("TARA_metadata/custom_tara_geo_metadata.txt", 
                        sep = "\t", header = T, as.is = T, fill = T, comment.char = "")
  colnames(geo_data)[15] = "Biogeographical_region"
  
  #Merge snp matrix with metadata
  meta_snps = merge(flipped_snps, geo_data, by.x = 0, by.y = "run_accession", all.x=T)
  
  #Calculate Principal Components
  snps.pca = prcomp(meta_snps[1:ncol(flipped_snps)])
  
  #Plot PCA
  autoplot(snps.pca, data=meta_snps, colour = "Biogeographical_region", legend = FALSE) +
    theme(legend.position = "none") + labs(title = pasate(s, " PCA, biogeographical regions colored", sep = ""))
  
}

# Close graphing device
dev.off()


















# setwd("/Users/nikaz/Desktop")
# 
# #Import snp frequency file (should be species specific)
# #gene_presence = read.table(file = "genes_presabs.txt", as.is = T, header = T)
# snp_freqs = read.table(file = "SFS_plots/snps_freqs/snps_ref_freq_macleodii.txt", as.is = T, header = T)
# 
# #Transpose data
# #flipped_genes = data.frame(t(gene_presence[-1]))
# #colnames(flipped_genes) = gene_presence[,1]
# 
# flipped_snps = data.frame(t(snp_freqs[-1]))
# colnames(flipped_snps) = snp_freqs[,1]
# 
# #Import metadata
# geo_data = read.table("TARA_metadata/custom_tara_geo_metadata.txt", 
#                       sep = "\t", header = T, as.is = T, fill = T, comment.char = "")
# colnames(geo_data)[15] = "Biogeographical_region"
# 
# sample = flipped_snps[1:3]
# 
# meta_snps = merge(flipped_snps, geo_data, by.x = 0, by.y = "run_accession", all.x=T)
# 
# library(ggfortify)
# #genes.pca = prcomp(flipped_genes)
# snps.pca = prcomp(meta_snps[1:ncol(flipped_snps)])
# 
# pdf("/Users/nikaz/Desktop/snp_PCA.pdf")
# 
# autoplot(snps.pca, data=meta_snps, colour = "Biogeographical_region", legend = FALSE) +
#   theme(legend.position = "none")
# 
# dev.off()
