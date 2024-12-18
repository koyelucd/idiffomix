#' @title Gene expression data for patients suffering from breast cancer
#'
#' @description A dataset containing simulated RNA-Seq data for \eqn{G} genes located on chromosomes 1 and 2, from \eqn{R=2} sample types, from \eqn{N=4} patients. The sample types are assumed to be benign and tumour tissues.
#' @details
#' The simulated raw RNA-Seq data for genes located on the chromosomes 1 and 2 needs to be filtered, normalized and transformed before applying idiffomix.
#'
#' @seealso \code{\link{gene_chromosome_data}}
#' @seealso \code{\link{methylation_data}}
#' @format A data frame with 20 rows and 9 columns. The data contain no missing values.
#' \itemize{
#'   \item{Gene: The gene name.}
#'   \item{Patient1_GX1: Expression values from benign tissue from patient 1.}
#'   \item{Patient2_GX1: Expression values from benign tissue from patient 2.}
#'   \item{Patient3_GX1: Expression values from benign tissue from patient 3.}
#'   \item{Patient4_GX1: Expression values from benign tissue from patient 4.}
#'   \item{Patient1_GX2: Expression values from tumour tissue from patient 1.}
#'   \item{Patient2_GX2: Expression values from tumour tissue from patient 2.}
#'   \item{Patient3_GX2: Expression values from tumour tissue from patient 3.}
#'   \item{Patient4_GX2: Expression values from tumour tissue from patient 4.}
#'    }
#' @usage data(gene_expression_data)
"gene_expression_data"
