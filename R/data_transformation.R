#' @title The function to filter, normalize and transform RNA-Seq and methylation data.
#' @export
#' @description The raw RNA-Seq and methylation data needs to be filtered, normalized and transformed before applying the idiffomix method.
#'
#' @details The RNA-Seq data consisted of raw counts depicting the gene expression levels. To ensure data quality, only genes whose sum of expression counts across both biological conditions was > 5 are retained. The data were normalized to account for differences in library sizes. The normalized count data were used to obtain CPM values which were further log-transformed to obtain log-CPM values. Given the paired design of the motivating setting, the log-fold changes between the tumour and benign samples were calculated for each gene in every patient and used in the subsequent analyses.
#' For the methylation array data, the beta values at the CpG sites are logit transformed to M-values. Similar to the RNA-Seq data, given the paired design, the difference in M-values between tumour and benign samples were calculated for each CpG site in every patient and used in the subsequent analyses.
#'
#' @param seq_data A dataframe of dimension \eqn{G \times (N+1)} containing raw RNA-Seq data for all G genes and N patients
#' @param meth_data A dataframe of dimension \eqn{C \times (N+2)} containing beta methylation values for all $C$ CpG sites and $N$ patients along with the associated genes for each CpG site.
#' @param gene_chr A dataframe containing the genes and their corresponding chromosome number.
#' @param N Number of patients
#'
#' @return The function returns a list with two dataframes containing the transformed gene expression and methylation array data:
#' \itemize{
#' \item seq_transformed - A dataframe containing the log-fold change for gene expression data.
#' \item meth_transformed - A dataframe containing the differences in M-values for methylation data.
#' }
#' @examples
#' N <- 4
#' data_output = data_transformation(seq_data=gene_expression_data,
#'                                   meth_data=methylation_data,
#'                                   gene_chr=gene_chromosome_data,
#'                                   N=N)
#'
#' @importFrom stats model.matrix
#' @importFrom edgeR filterByExpr DGEList calcNormFactors cpm
#' @importFrom magrittr %>%
#'

data_transformation<-function(seq_data,meth_data,gene_chr,N=5){

    ### Filtering, normalizing and transforming gene expression data

    rownames(seq_data)<-seq_data$Gene
    gene_df_new<-seq_data[,-1]
    targets <- data.frame(Patients=c(seq(1,N),seq(1,N)),
                          Conditions=c(rep("Benign", N),rep("Tumour", N)))
    Patients <- factor(targets$Patients)
    Conditions <- factor(targets$Conditions, levels=c("Benign","Tumour"))
    design <- model.matrix(~0+Patients+Conditions)
    design

    keep <- filterByExpr(gene_df_new, design)
    gene_df_new2 <- gene_df_new[keep,]
    gene_cpm<-gene_df_new2
    # Create DGEList object
    dge <- DGEList(counts = gene_cpm)
    nonzero <- rowSums(dge$counts) > 5
    # dge %>% .[nonzero,]
    dge[nonzero,]
    dge %>% calcNormFactors
    cpm_values <- cpm(dge)

    # Transform filtered RNA-Seq count data

    log_cpm_values <- log2(cpm_values + 1)  # Adding 1 to avoid log of zero
    log_cpm_values<-as.data.frame(log_cpm_values)
    log_cpm_values<-log_cpm_values[order(rownames(log_cpm_values)),]
    seq_transformed_df<-log_cpm_values
    seq_transformed_df<-as.data.frame(cbind(rownames(log_cpm_values),seq_transformed_df))
    colnames(seq_transformed_df)[1]<-"Gene"
    cols_to_convert <- colnames(seq_transformed_df)[-1]
    seq_transformed_df[cols_to_convert] <- lapply(seq_transformed_df[cols_to_convert], as.numeric)

    for (i in 2:(N+1)) {
        seq_transformed_df[, i] <- seq_transformed_df[, i + N] - seq_transformed_df[, i]
    }
    seq_transformed_df <- seq_transformed_df[, -((N+2):(2*N+1))]
    seq_transformed_df<-as.data.frame(seq_transformed_df)
    colnames(seq_transformed_df) <- c("Gene", paste0("Patient", 1:N))
    cols_to_convert <- paste0("Patient", 1:N)
    seq_transformed_df[cols_to_convert] <- lapply(seq_transformed_df[cols_to_convert], as.numeric)



    #### Transforming methylation array data
    meth_transformed_df<-meth_data
    meth_transformed_df[,(3:((N*2)+2))]<-log(meth_transformed_df[,(3:((N*2)+2))]/(1-meth_transformed_df[,(3:((N*2)+2))]))
    for (i in 3:(N+2)) {
        meth_transformed_df[, i] <- meth_transformed_df[, i + N] - meth_transformed_df[, i]
    }
    meth_transformed_df <- meth_transformed_df[, -((N+3):(2*N+2))]

    colnames(meth_transformed_df)<-c("Gene","CpG",paste0("Patient", 1:N))
    cols_to_convert <- paste0("Patient", 1:N)
    meth_transformed_df[cols_to_convert] <- lapply(meth_transformed_df[cols_to_convert], as.numeric)

    seq_transformed_df<-seq_transformed_df[seq_transformed_df$Gene %in% gene_chr$Gene,]
    meth_transformed_df<-meth_transformed_df[meth_transformed_df$Gene %in% gene_chr$Gene,]

    return(list(seq_transformed=seq_transformed_df,
                meth_transformed=meth_transformed_df))

}
