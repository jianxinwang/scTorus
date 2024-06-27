#' Data from Vectra for a cell cycle panel 
#'
#' Shrimp are classified by size, 0-15 shrimp per pound, 15-20 shrimp per pound, etc. A smaller number per pound indicates larger shrimp. Nominal prices are total monthly value of brown shrimp andings within size class divided by total monthly landings within the size class. 
#'
#' @format A data table with 275709 rows and 10 columns:
#' \describe{
#'   \item{`Sample Name`}{char TMA block and ROI IDs}
#'   \item{`Cell ID`}{int Numberic Cell ID for each ROI} 
#'   \item{`Nucleus DAPI Mean`}{dbl denoting scanned intensity for nucleus DAPI }
#'   \item{`Nucleus CCNE Mean`}{dbl denoting scanned intensity for nucleus Cyclin E1}
#'   \item{`...`}{dbl Other markers...}
#' }
#' 
#' @source Test dataset from lab
#'
#' @examples
#' data(cell_cycle)
"cell_cycle"