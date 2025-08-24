
library(tidyverse)


region <- "22_24234172_24401503" # MIF
region <- "22_23942068_26128449" # GGT1

path_base <- "/scratch/dariush.ghasemi/projects/susie_pqtls/ld_block/results/"
path_ld <- paste0(path_base, "ld/", region, ".ld")
path_ld_matrix <- paste0(path_base, "ldmat/", region, ".ld")
path_ld_snps <- paste0(path_base, "ldmat/", region, ".snps")

# Load LD matrix
# change class to matrix
ld_matrix <- as.matrix(data.table::fread(path_ld_matrix))

# Load SNP positions
snp_positions <- read.table(path_ld_snps, col.names = c("SNP", "Position"))

# Prepare heatmap
colnames(ld_matrix) <- snp_positions$Position
rownames(ld_matrix) <- snp_positions$Position


pdf("LD_blocks_GGT1.pdf", width = 15, height = 15) #, units = "cm", res = 300
jpeg("LD_blocks_GGT1.jpg", width = 15, height = 15, units = "cm", res = 700)

# Simple heatmap with base R:
heatmap(
  ld_matrix, 
  Rowv = NA, 
  Colv = NA, 
  scale = "none",
  main = "GGT1 region",
  col = colorRampPalette(c("white", "red"))(100)
  )

# add legend to heatmap
# legend(
#   x= 0.9, y = 0.5, #"bottomright", 
#   legend=c("0", "0.5", "1"),
#   fill= colorRampPalette(c("white", "red"))(3)
#   )


dev.off()




ld_long <- data.table::fread(path_ld) #%>%
#  as.data.frame() %>%
#  rownames_to_column("SNP1") %>%
#  pivot_longer(-SNP1, names_to = "SNP2", values_to = "r2")

# Or better: ggplot2 heatmap
ld_long %>%
  dplyr::mutate(across(c("BP_A", "BP_B"), as.character)) %>% 
  ggplot(aes(x = BP_A, y = BP_B, fill = R2)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "LD Heatmap", x = "SNP", y = "SNP")



#=============================#
# -----    MIF index     -----
#=============================#

library(corrplot)

# calculate LD between index variants of the 6 loci in broader MIF region 
mif_index <- "/scratch/dariush.ghasemi/projects/susie_pqtls/ld_block/results/dosage/mif_cojo.raw"

#---- Read genotype file ----
genos <- fread(mif_index)

# remove effective allele from column names
genos <- genos %>% rename_with(~ str_remove(., "_[ATCG]+$"), matches("_[ATCG]+$"))

# store variants list for the match
snps <- names(genos)[7:ncol(genos)]

X <- as.matrix(genos[,7:ncol(genos)]) * 1.0

# compute LD (r2)
R <- cor(X, use = "complete.obs")  # ---> Find a solution for NAs in genotype and cor matrix


corrplot.mixed(R)

corrplot(
  R, 
  method = 'number', 
  diag = FALSE,
  #type = 'upper',
  col = COL2('PRGn') #COL2('BrBG')
  ) # colorful number



# Gene is used for annotation of Pheatmap
annot_heatmap <- data.frame(
  row.names = colnames(R),
  Gene = c("MIF", "MIF", "DDT", "GSTT1", "GSTT1", "DDT", "GSTT1", "GGT5", "GGT5", "GGT5", "F10")
)

pheatmap(R, annotation_col = annot_heatmap)

