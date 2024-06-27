## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

library(data.table)
library(dplyr)

# BCM data
setwd("/mnt/Vectra/Vectra Projects/Vectra Projects_ Witkiewicz.Knudsen Lab/2023_9 EK Baylor TMA/Consolidated for Jason/")

data.bcm.cc2.09.1 <- fread("Consolidated_data_BCM09-1.txt")
data.bcm.cc2.09.2 <- fread("Consolidated_data_BCM09-2.txt")
data.bcm.cc2.10.1 <- fread("Consolidated_data_BCM10-1.txt")
data.bcm.cc2.10.2 <- fread("Consolidated_data_BCM10-2.txt")

tmpList <- list()
tmpList[['bcm.cc2.09.1']] <- data.bcm.cc2.09.1
tmpList[['bcm.cc2.09.2']] <- data.bcm.cc2.09.2
tmpList[['bcm.cc2.10.1']] <- data.bcm.cc2.10.1
tmpList[['bcm.cc2.10.2']] <- data.bcm.cc2.10.2

dataSetList <- list()
getwd()
setwd('/projects/Rpackages/scTorus/inst/extdata')
for (ds in names(tmpList)){
    tmp <- tmpList[[ds]]
    tmp <- tmp[`Tissue Category` == 'Tumor']
    
    # get the mean value columns
    tmp <- tmp %>% dplyr::select(contains(c('Sample Name', 'Cell ID', 'Mean'))) %>% 
        dplyr::select(-contains(c('AE1AE3', 'Autofluorescence', 'Nucleus', 'Membrane', 'Cytoplasm')))

    
    # clean column names
    colnames(tmp) <- sub(' \\(Opal \\d+\\)', '', colnames(tmp))
    dataSetList[[ds]] <- tmp
    
    # save to file
    saveRDS(tmp, file=paste0(ds, '.rds'))
}





