## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = TRUE)
#knitr::opts_chunk$set(dev = 'png')
#knitr::opts_chunk$set(dpi=100)

## ----biocsetup, eval=FALSE----------------------------------------------------
#  library(plotly)
#  library(ggplot2)
#  library(foreach)
#  library(doParallel)
#  
#   # library(devtools)
#   # install_github('koyelucd/idiffomix',force = TRUE)
#   # library(idiffomix)

## ----package, include=TRUE, echo=TRUE, message=FALSE,warning=FALSE------------
library(idiffomix)

## ----data,include=TRUE, echo=TRUE---------------------------------------------
data(gene_expression_data)
data(methylation_data)
data(gene_chromosome_data)

## ----data_transformation, include=TRUE, echo=TRUE-----------------------------
N <- 4
data_transformed = data_transformation(seq_data=gene_expression_data,
                                       meth_data=methylation_data,
                                       gene_chr=gene_chromosome_data,
                                       N=N)

## ----idiffomix, include=TRUE, echo=TRUE---------------------------------------
N <- 4
data_output = idiffomix(seq_data=data_transformed$seq_transformed,
                        meth_data=data_transformed$meth_transformed,
                        gene_chr=gene_chromosome_data,
                        N=N, K=3, L=3, probs=c(0.25,0.75),
                        parallel_process = FALSE)

print(data_output$tau[[1]])
head(data_output$seq_classification)

## ----chromosome_plot, include=TRUE, echo=TRUE,fig.width=8, fig.height=6,dev = 'png',warning=FALSE----
plot(data_output,what="chromosome",CHR=1, Gene=NULL,K=3,L=3,
     gene_cluster_name=c( "E-","E0","E+"),
     cpg_cluster_name=c( "M-","M0","M+"))

## ----gene_plot, include=TRUE, echo=TRUE,fig.width=14, fig.height=8,dev = 'png',warning=FALSE----
plot(data_output,what="gene",CHR=1, Gene="gene_1",K=3,L=3,
     gene_cluster_name=c( "E-","E0","E+"),
     cpg_cluster_name=c( "M-","M0","M+"))

