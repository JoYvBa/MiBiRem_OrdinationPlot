### Fill in path to the Sample_data directory to make it the working directory ###
setwd("~/Sample_data")
##################################################################################
Spec_sheet <- read.csv("CCA_validation_species.csv", row.names=1)
Env_sheet <- read.csv("CCA_validation_environment.csv", row.names=1)
cca_res <- cca(Spec_sheet, Env_sheet, scale=FALSE)
species <- cca_res[["CCA"]][["v"]]
enviroment <- cca_res[["CCA"]][["biplot"]]
species_12 <- as.data.frame(species[,1:2])
enviroment_12 <- as.data.frame(enviroment[,1:2])
write_xlsx(species_12, "CCA_validation_species_loadings.xlsx")
write_xlsx(enviroment_12, "CCA_validation_environment_loadings.xlsx")
