setwd("/Users/nikaz/Downloads")

# Read in file containing the top 30 lines of each species profile (from every TARA sample)
# species = read.table(file="top_30_species.txt", as.is = T)
species = read.table(file="/Users/nikaz/Desktop/species_merge.txt", as.is = T, header = T)
# colnames(species) = c("species_id", "count_reads", "coverage", "relative_abundance") # Assign column names

species = subset(species, count_reads != 0) # Remove all entries with zero coverage. 

head_1 = read.table(file="/Users/nikaz/Desktop/top_species.txt", as.is = T, header = T)

#This command takes just the species ID column from the original dataset, makes a table of counts
#(how many times each species name appears), then sorts the table by count in decreasing order.
species_counts = sort(table(species$species_id), decreasing = T) 

#Plotting
colors = rep("tomato", ) #This attempt at coloring the subsequent barplot is not working.
barplot(species_counts, horiz=T, col = colors) #Make a barplot with species_counts

#Saving results
species_table = as.data.frame(species_counts) 
colnames(species_table) = c("Species", "Count")

write.csv(species_table[1:20,], file = "top_20_species.csv", row.names = F)



## Coverage as a function of Counts ##
group = factor(species$species_id)
plot(species$count_reads, species$coverage, 
     main = "Read count vs Coverage", xlab = "Reads mapped to target genes", ylab = "Read Depth",
     col = group, pch = 16)
## Legend code, (its to big to be practical) ##
#legend("bottomright", 
#       legend = levels(group), col = unique(group)


##Coverage Distribution##

top_20 = species_table[1:20,]
top_20_coverage = subset(species, species_id %in% top_20$Species, c(species_id, coverage))
hist(top_20_coverage$coverage, plot=T, main = "Top 20 Species Coverage Distribution",
     col = "peachpuff", xlab = "coverage", freq = T, breaks=100)
rug(top_20_coverage$coverage)

hist(species$coverage, plot=T, main = "All Samples Coverage Distribution",
     col = "green", xlab = "coverage", freq = T, breaks=100)
rug(top_20_coverage$coverage)

hist(head_1$coverage, plot=T, main = "Top per Sample Species Coverage",
     col = "green", xlab = "coverage", freq = T, breaks=20)
rug(head_1$coverage)
abline(v=142, pch=500)


bplot = boxplot(head_10$coverage)
outliers = which(top_20_coverage$coverage %in% bplot$out)

top_20_coverage[outliers,]
