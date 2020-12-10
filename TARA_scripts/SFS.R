setwd("/Users/nikaz/Desktop")

#Import snp frequency file (should be species specific)
snp_freqs = read.table(file = "SFS_plots/snps_freqs/snps_ref_freq_macleodii.txt", as.is = T, header = T)

#Delete unnecessary column
snp_freqs = subset(snp_freqs, select = -(site_id))

## SFS Plot per Sample ##

prep_snps = function(snps){
  #Filter snp frequencies by these min/max thresholds
  trimmed = snps[snps < 0.95 & snps > 0.05] 
  #Fold the plot (any frequencies below 0.5 get subtracted from 1)
  for (freq in 1:length(trimmed)){
    if (trimmed[freq] < 0.5){
      trimmed[freq] = 1 - trimmed[freq]
    }
  }
  return(trimmed)
}

#Open plotting device (pdf)
pdf("/Users/nikaz/Desktop/snp_freqs.pdf")

#Loop for plotting one SFS per sample
for (sample in 1:ncol(snp_freqs)){
  #Extract one sample's worth of snps
  snps = snp_freqs[,sample]
  #Filter snp frequencies by these min/max thresholds
  trimmed = snps[snps < 0.95 & snps > 0.05] 
  #Fold the plot (any frequencies below 0.5 get subtracted from 1)
  for (freq in 1:length(trimmed)){
    if (trimmed[freq] < 0.5){
      trimmed[freq] = 1 - trimmed[freq]
    }
  }
  #Extract sample name
  sample_name = colnames(snp_freqs)[sample]
  #Plot histogram
  hist(trimmed, freq = T, breaks = 30, col="magenta", 
       main = paste("B_vulgatus Sample ", sample_name, " SFS"), 
       xlab = "Allele Frequency", xlim = c(0.45, 1))
} 

#Close graphing device
dev.off()


## Plot ranked samples based on mid-range frequencies ##

#Calculate the fraction of snps with frequencies between 0.2 and 0.8
mid_freqs = apply(snp_freqs, 2, function(x) sum(x <= 0.8 & x >= 0.2)/length(x))
#Sort in descending order
mid_freqs = sort(mid_freqs, decreasing = T)

#Code for coloring specific samples
colors = rep("black", length(mid_freqs))
#colors = rainbow(length(mid_freqs))

for (i in 1:length(colors)){
  if (names(mid_freqs)[i] == "ERR315857"){
    colors[i] = "black"
  }else if (names(mid_freqs)[i] == "ERR599043"){
    colors[i] = "black"
  }else if (names(mid_freqs)[i] == "ERR594345"){
    colors[i] = "red"
  }
}

#Plot scatter plot
plot(mid_freqs, col = colors, pch = 19, 
     main = "A_macleodii",
     ylab = "Within sample polymorphism (0.2 ≤ f ≤ 0.8)",
     xlab = paste("Ranked Samples, n = ", length(mid_freqs))

)



  

## SFS Plot per Geographic Variable ##
#Import TARA metadata (see TARA_metadata.R for how this is organized)
# NOTE: the header contains a "#" symbol which R reads as a commnet -> use comment.char = ""
geo_data = read.table("TARA_metadata/custom_tara_geo_metadata.txt", 
                      sep = "\t", header = T, as.is = T, fill = T, comment.char = "")

#Combine samples by geographical province (column 15 in geo_data)
#stack function useful to combine columns, stack(snp_freqs, select = samples in region)

regions = levels(factor(geo_data[,15]))
colnames(geo_data)[15] = "Biogeographical_region"

#Open plotting device (pdf)
pdf("/Users/nikaz/Desktop/snp_freqs.pdf")

for (region in regions){
  #Extract sample accessions which exist in this species and were taken from specific region
  samples = subset(geo_data, Biogeographical_region == region & run_accession %in% colnames(snp_freqs), 
                   select = run_accession)
  #If any samples fit the criteria, stack their snps into one column, and plot
  if(length(samples$run_accession)){
    snps = stack(snp_freqs, select = c(samples$run_accession))
    #Change snp freqs to a vector
    snps = snps$values
    #Prep snps for plotting
    trimmed = prep_snps(snps)
    
    #Plot histogram
    hist(trimmed, freq = T, breaks = 30, col="tomato", 
         main = paste("A_macleodii in\n", region, "\n", length(samples$run_accession), "combined samples"), 
         xlab = "Allele Frequency", xlim = c(0.45, 1))
    

  }
}
#Close plotting device
dev.off()

# Combine the two most popolous regions (Indian Monssoon and Mediterranean Sea)
regions = levels(factor(geo_data[,15]))
most_prevalent_regions = regions[c(12, 13)]
samples = subset(geo_data, Biogeographical_region %in% most_prevalent_regions & run_accession %in% colnames(snp_freqs), 
                 select = run_accession)
snps = stack(snp_freqs, select = c(samples$run_accession))
snps = snps$values
trimmed = prep_snps(snps)
hist(trimmed, freq = T, breaks = 30, col="orange", 
     main = paste("A_macleodii in most prevalent regions\n", length(samples$run_accession), "combined samples"), 
     xlab = "Allele Frequency", xlim = c(0.45, 1))
