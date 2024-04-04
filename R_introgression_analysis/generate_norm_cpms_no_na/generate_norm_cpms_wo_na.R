library(changepoint)

load("/home/adejoro/Desktop/R_folders/change_point_2/norm_lyc_cpm")
norm_lyc_cpm <- cpm_filtered
rm(cpm_filtered)

load("/home/adejoro/Desktop/R_folders/change_point_2/lyc_cpm")
lyc_cpm <- cpm_filtered
rm(cpm_filtered)
any(is.nan(lyc_cpm))

###
# Calculate row means
row_means <- rowMeans(lyc_cpm)

# Replace row means that are zero with a small positive value (e.g., 1e-10)
row_means[row_means == 0] <- 1e-10

# Divide each element of cpm_filtered by the corresponding row mean
norm_lyc_cpm_no_nas <- lyc_cpm / row_means
save(norm_lyc_cpm_no_nas, file="norm_lyc_cpms_wo_na")

###
pi <- which(colnames(norm_lyc_cpm_no_nas) == "SRR11874184")
cur_samp <- norm_lyc_cpm_no_nas[, pi]

plot(log(cur_samp))

###
rm(lyc_cpm)
rm(norm_lyc_cpm)

###
load("/home/adejoro/Desktop/R_folders/change_point_4/pen_cpm")
pen_cpm <- cpm_filtered
rm(cpm_filtered)

row_means_pen <- rowMeans(pen_cpm)
row_means_pen[row_means_pen == 0] <- 1e-10
norm_pen_cpms_no_nas <- pen_cpm / row_means_pen
save(norm_pen_cpms_no_nas, file="norm_pen_cpms_wo_na")

pi <- which(colnames(norm_pen_cpms_no_nas) == "SRR11874184")
cur_samp <- norm_pen_cpms_no_nas[, pi]

plot(log(cur_samp))
