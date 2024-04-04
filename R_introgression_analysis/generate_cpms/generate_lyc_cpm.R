library(tidyverse)
library(edgeR)

lyc_matrix_filtered <- read.csv("/home/adejoro/Desktop/from-server/filtered_counts_matrix/combined_filtered_counts_matrix/combined_filtered_matrix.csv")
save(lyc_matrix_filtered, file = "untransposed_lyc_filtered_matrix")

lyc_matrix_filtered <- t(lyc_matrix_filtered)

colnames(lyc_matrix_filtered) <- lyc_matrix_filtered[1,]
lyc_matrix_filtered <- lyc_matrix_filtered[-1, ]

filtered_tibble <- as_tibble(lyc_matrix_filtered, rownames="geneID")

filtered_geneIDS <- filtered_tibble$geneID
save(filtered_geneIDS, file="lyc_gene_IDs")

filtered_tibble <- apply(filtered_tibble[, -1], 2, as.numeric)
save(filtered_tibble, file="filtered_lyc_tibble")


lyc_matrix_filtered <- as.matrix(filtered_tibble)
save(lyc_matrix_filtered, file="lyc_matrix_filtered")

filtered_DGE <- DGEList(lyc_matrix_filtered)
save(filtered_DGE, file="filtered_lyc_DGE")

filtered_DGE <- calcNormFactors(filtered_DGE, method = "TMM")
save(filtered_DGE, file="norm_fil_lyc_DGE")

cpm_filtered <- cpm(filtered_DGE, log=FALSE)
save(cpm_filtered, file="lyc_cpm")
dim(cpm_filtered)

cpm_filtered <- cpm_filtered / rowMeans(cpm_filtered)
dim(cpm_filtered)
save(cpm_filtered, file="norm_lyc_cpm")
