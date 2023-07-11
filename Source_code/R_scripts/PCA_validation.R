### Fill in path to the Sample_data directory to make it the working directory ###
setwd("~/Sample_data")
##################################################################################
PCA_sheet <- read.csv("PCA_validation.csv", row.names=1)
pca_res <- prcomp(PCA_sheet, scale = FALSE)
loadings <- pca_res$rotation
loadings_12 <- as.data.frame(loadings[,1:2])
write_xlsx(loadings_12, "PCA_validation_loadings.xlsx")
