## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

library(data.table)
library(dplyr)

# BCM data

setwd("/mnt/Vectra/Vectra_Projects/Vectra_Projects_Witkiewicz.Knudsen_Lab/2023.9 EK Baylor TMA/Consolidated for Jason/")

data.bcm.cc2.09.1 <- fread("Consolidated_data_BCM09-1.txt")
data.bcm.cc2.10.2 <- fread("Consolidated_data_BCM10-2.txt")

tmpList <- list()
tmpList[['cell_cycle']] <- data.bcm.cc2.09.1


dataSetList <- list()
getwd()
setwd('/projects/Rpackages/scTorus/inst/extdata')
for (ds in names(tmpList)){
    tmp <- tmpList[[ds]]
    tmp <- tmp[`Tissue Category` == 'Tumor']
    
    # get the mean value columns
    tmp <- tmp %>% dplyr::select(contains(c('Sample Name', 'Cell ID', 'Mean'))) %>% 
        dplyr::select(-contains(c('AE1AE3', 'Autofluorescence', 'DAPI', 'Entire', 'Membrane', 'Cytoplasm')))

    
    # clean column names
    colnames(tmp) <- sub(' \\(Opal \\d+\\)', '', colnames(tmp))
    dataSetList[[ds]] <- tmp
    
    # save to file
    saveRDS(tmp, file=paste0(ds, '.rds'))
}


usethis::use_r("cell_cycle")

devtools::document()
devtools::check()
