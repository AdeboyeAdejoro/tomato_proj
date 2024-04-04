library(edgeR)
library(tidyverse)

pen_mat_filtered <- read.csv("/home/adejoro/Desktop/from-server/filtered_counts_matrix_pen/filtered_counts_pen_combined_clean/filtered_pen_counts.csv")

pen_mat_filtered <- t(pen_mat_filtered)

colnames(pen_mat_filtered) <- pen_mat_filtered[1,]
pen_mat_filtered <- pen_mat_filtered[-1, ]

filtered_tibble <- as_tibble(pen_mat_filtered, rownames="geneID")

filtered_geneIDS <- filtered_tibble$geneID
save(filtered_geneIDS, file="pen_geneIDs")

filtered_tibble <- apply(filtered_tibble[, -1], 2, as.numeric)

pen_matrix_filtered <- as.matrix(filtered_tibble)
save(pen_matrix_filtered, file="filtered_pen_matrix_v1")

filtered_DGE <- DGEList(pen_matrix_filtered)

filtered_DGE <- filtered_DGE[, filtered_DGE$samples$lib.size > 0]
save(filtered_DGE, file="filtered_filtered_DGE")

filtered_DGE <- calcNormFactors(filtered_DGE, method = "TMM")
save(filtered_DGE, file="norm_DGE")

cpm_filtered <- cpm(filtered_DGE, log=FALSE)
save(cpm_filtered, file="pen_cpm")

cpm_filtered <- cpm_filtered / rowMeans(cpm_filtered)
save(cpm_filtered, file="norm_cpm")


