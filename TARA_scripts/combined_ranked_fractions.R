
## Combined ranks from all species ##

setwd("/u/scratch/n/nzeltser/midas/snps")

# Read in file containing name of every species 
species = read.table(file = "species_names.txt", header = T, as.is = T)

# Initialize empty vector (to be filled with polymorphism fractions)
poly_fractions = vector()

# Iterate over each species
for (s in species$species_name){
  # read in snp frequency file
  snp_freqs = read.table(file = paste(s, "/snps_ref_freq.txt", sep = ""), header = T, as.is = T)
  snp_freqs = subset(snp_freqs, select = -(site_id))
  # read in snp coverage file
  snp_coverage = read.table(file = paste(s, "/snps_depth.txt", sep = ""), header = T, as.is = T)
  snp_coverage = subset(snp_coverage, select = -(site_id))
  # Remove snps with lower than 20x coverage from future calculations by changing their frequency to NA
  # This expression works because the order of snps and samples should be identical in the two files (so indexing matches up)
  snp_freqs[snp_coverage < 20] = NA
  # For each sample, add up the number of snps whose frequency falls within the threshold, then divide by the total number of snps (excluding NA values)
  mid_freq = apply(snp_freqs, 2, function(x) sum(x <= 0.8 & x >= 0.2, na.rm = T)/sum(!is.na(x)))
  # Append new fractions to existing vector
  poly_fractions = c(poly_fractions, mid_freq)
}

# Sort fractions by decreasing order
poly_fractions = sort(poly_fractions, decreasing = T)

# Open link to pdf plotting device
pdf("ranked_polymorphism.pdf")
# Plot as a ranked scatter plot, with log transformed y-axis
plot(poly_fractions, col = "black", pch = 19, log = "y",
     ylab = "Within sample polymorphism (0.2 ≤ f ≤ 0.8)",
     xlab = paste("Ranked Samples, n = ", length(poly_fractions))
     
)

# Close plotting device
dev.off()






# ## Combined ranks from all species ##
# 
# species = read.table(file = "species_name.txt", header = T, as.is=T)
# 
# #species = c("snps_ref_freq_macleodii.txt", "snps_ref_freq_marinus.txt")
# poly_fractions = vector()
# 
# for (s in species$species_name){
#   snp_freqs = read.table(file = paste("SFS_plots/snps_freqs/snps_ref_freq_",s,".txt", sep = ""), header = T, as.is = T)
#   snp_freqs = subset(snp_freqs, select = -(site_id))
#   snp_coverage = read.table(file = paste("SFS_plots/snps_depth/snps_depth_",s,".txt", sep = ""), header = T, as.is = T)
#   snp_coverage = subset(snp_coverage, select = -(site_id))
#   snp_freqs[snp_coverage < 20] = NA
#   mid_freq = apply(snp_freqs, 2, function(x) sum(x <= 0.8 & x >= 0.2, na.rm = T)/sum(!is.na(x)))
#   poly_fractions = c(poly_fractions, mid_freq)
# }
# 
# poly_fractions = sort(poly_fractions, decreasing = T)
# 
# pdf("ranked_polymorphism.pdf")
# plot(poly_fractions, col = "black", pch = 19,
#      ylab = "Within sample polymorphism (0.2 ≤ f ≤ 0.8)",
#      xlab = paste("Ranked Samples, n = ", length(poly_fractions))
#      
# )
# dev.off()
