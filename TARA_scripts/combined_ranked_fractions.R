## Combined ranks from all species ##

species = read.table(file = "species_name.txt", header = T)

#species = c("snps_ref_freq_macleodii.txt", "snps_ref_freq_marinus.txt")
poly_fractions = vector()

for (s in species$species_name){
  snp_freqs = read.table(file = paste(s, "snps_ref_freq.txt", sep = ""), header = T, as.is = T)
  snp_freqs = subset(snp_freqs, select = -(site_id))
  mid_freq = apply(snp_freqs, 2, function(x) sum(x <= 0.8 & x >= 0.2)/length(x))
  poly_fractions = c(poly_fractions, mid_freq)
}

poly_fractions = sort(poly_fractions, decreasing = T)

pdf("ranked_polymorphism.pdf")
plot(poly_fractions, col = "black", pch = 19,
     ylab = "Within sample polymorphism (0.2 ≤ f ≤ 0.8)",
     xlab = paste("Ranked Samples, n = ", length(poly_fractions))
     
)
dev.off()


