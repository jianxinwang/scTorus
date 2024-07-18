#' make_input_from_vectra
#'
#' This function takes a consolidated data file from Vectra and output a data frame
#' as input for the function of run_scTorus for cell cycle trajectory analysis
#' and trajectory profiles showing the up- and downs of various cell cycle protine markers.
#' @param infile A tab delimited text file from Vectra machine
#' @param tissue_category What type of cells we are interested in. Options are "Tumor", "Stroma"
#' @param compartment Intensity values measured from which cellular compartment. Options are
#' 'Nuleus', 'Entire.Cell' and 'Cytoplasm'
#' @param excluded_markers Marker(s) to exlude from the output. For example, we may not want 'AE1AE3' and/or 'DAPI'
#' @return A data frame for use as input for the run_scTorus function
#' @export
#'

make_input_from_vectra <- function(infile, tissue_category = c('Tumor', 'Stroma'),
                                   compartment = c('Nucleus', 'Entire.Cell', 'Cytoplasm'),
                                   excluded_markers = c('AE1AE3', 'DAPI')){
    
    dat <- data.table::fread(infile)
    dat <- as.data.frame(dat)
    
    # Fix column naming inconsistencies
    colnames(dat) <- gsub("[ |\\(|\\)]", '.', colnames(dat))
    
    rownames(dat) <- paste(dat$Sample.Name, dat$Cell.ID, sep = ':')
    
    # Only work on tumor cells
    dat.tumor <- dat %>% filter(Tissue.Category == tissue_category)
  
    if (compartment == 'Nucleus'){
        data.for.trajectory <- dat.tumor %>% dplyr::select(matches('Nucleus.*Mean'))
        colnames(data.for.trajectory) <- sub('^Nucleus\\.(.*?)\\..*', '\\1', colnames(data.for.trajectory))
    } else if (compartment == 'Entire.Cell'){
        data.for.trajectory <- dat.tumor %>% dplyr::select(matches('Entire.*Mean'))
        colnames(data.for.trajectory) <- sub('^Entire\\.Cell\\.(.*?)\\..*', '\\1', colnames(data.for.trajectory))
    } else if (compartment == 'Cytoplasm'){
        data.for.trajectory <- dat.tumor %>% dplyr::select(matches('Cytoplasm.*Mean'))
        colnames(data.for.trajectory) <- sub('^Cytoplasm\\.(.*?)\\..*', '\\1', colnames(data.for.trajectory))
    } else {
        stop(paste("Unknown compartment:", compartment))
    }
    
    colnames(data.for.trajectory) <- sub('\\+', '', colnames(data.for.trajectory) )
    
    
    # Drop the unwanted columns
    data.for.trajectory <- data.for.trajectory %>% dplyr::select(-any_of(c('AF', 'Autofluorescence', excluded_markers)))
    
    data.for.trajectory
}

