### Fill in path to the Sample_data directory to make it R working directory ###
setwd("~/Sample_data")
##################################################################################
Spec_sheet <- read.csv("RDA_validation_species.csv", row.names=1)
Env_sheet <- read.csv("RDA_validation_environment.csv", row.names=1)
rda_res <- rda(Spec_sheet, Env_sheet, scale=FALSE)
species <- rda_res[["CCA"]][["v"]]
enviroment <- rda_res[["CCA"]][["biplot"]]
species_12 <- as.data.frame(species[,1:2])
enviroment_12 <- as.data.frame(enviroment[,1:2])
write_xlsx(species_12, "RDA_validation_species_loadings.xlsx")
write_xlsx(enviroment_12, "RDA_validation_environment_loadings.xlsx")
