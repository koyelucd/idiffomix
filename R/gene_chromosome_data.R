#' @title Data containing chromosome information and the genes located on them.
#'
#' @description A dataset containing the chromosome information of the gene expression and methylation array data to be analysed.
#' @details
#' The dataset contains the information of chromosomes 1 AND 2 and the genes located on them.
#'
#' @seealso \code{\link{gene_expression_data}}
#' @seealso \code{\link{methylation_data}}
#' @format A data frame with 20 rows and 2 columns.
#' \itemize{
#'   \item{CHR: The chromosome containing the gene.}
#'   \item{Gene: The gene located on the chromosome.}
#'    }
#' @usage data(gene_chromosome_data)
"gene_chromosome_data"
