utils::globalVariables(c(
    "Cluster", "CpG", "CpG_Cluster", "Final_Cluster", "Gene_Cluster",
    "Probability", "Type", "Value", "text_color"
))
#' @title Plots for visualizing the idiffomix class object
#' @description Visualise a \code{\link[idiffomix:idiffomix]{idiffomix}} clustering solution by plotting the conditional probabilities estimated for genes and CpG sites (A) per chromosome and (B) for a gene and its corresponding CpG sites.
#' @details The function displays two plots. The first plot displays the conditional probabilities estimated when the joint mixture model is applied to data from a chromosome.
#' Panel A displays the probability of a gene in the chromosome belonging to each of the K clusters.
#' Panel B details the estimated matrix pi of conditional probabilities of a CpG site's cluster membership, given its gene's cluster.
#' Panel C details the conditional probabilities of a gene belonging to cluster $k$ given a single CpG site associated with the gene belongs to cluster $l$, computed using Bayes' theorem,
#' given tau and pi. The second plot displays the log-fold changes and differences in M-values and the estimated posterior probability of the gene belonging to the K clusters.
#' Panel A shows the log-fold change and difference in M-values for the given gene and its associated CpG sites while Panel B shows the posterior probabilities of cluster membership for the gene under idiffomix.
#' @rdname plot.idiffomix
#' @export
#'
#' @seealso \code{\link{idiffomix}}
#' @param x A \code{\link[idiffomix:idiffomix]{idiffomix}} object.
#' @param what The different plots that can be obtained are either "chromosome" or "gene" (default = "chromosome").
#' @param CHR The chromosome number to be visualized (default = 1).
#' @param Gene The name of the gene to be visualized (default = NULL).
#' @param K Number of clusters for expression data (default = 3; E-,E0,E+).
#' @param L Number of clusters for methylation data (default = 3; M-,M0,M+).
#' @param gene_cluster_name The names to be given to the clusters for identification (default = c( "E-","E0","E+")).
#' @param cpg_cluster_name The names to be given to the clusters for identification (default = c( "M-","M0","M+")).
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user:
#' \itemize{
#' \item chromosome plot - Plot showing the conditional probabilities estimated when the joint mixture model is applied to data from a chromosome.
#' \item gene plot - Plot showing the log-fold changes and differences in M-values and the estimated posterior probability of the gene belonging to the K clusters.
#' }
#'
#' @examples
#' N <- 4
#' data_transformed = data_transformation(seq_data=gene_expression_data,
#'                                   meth_data=methylation_data,
#'                                   gene_chr=gene_chromosome_data,
#'                                   N=N)
#' data_output = idiffomix(seq_data=data_transformed$seq_transformed,
#'                         meth_data=data_transformed$meth_transformed,
#'                         gene_chr=gene_chromosome_data,
#'                         N=N, K=3, L=3, probs=c(0.25,0.75),
#'                         parallel_process = FALSE)
#' plot(data_output,what="chromosome",CHR=1, Gene=NULL,K=3,L=3,
#'      gene_cluster_name=c( "E-","E0","E+"),
#'      cpg_cluster_name=c( "M-","M0","M+"),
#'      title=NULL)
#'
#' @importFrom ggplot2 ggplot geom_bar scale_fill_gradient xlab ylim aes ylab
#' @importFrom ggplot2 theme_minimal theme element_text geom_text geom_tile
#' @importFrom ggplot2 geom_boxplot geom_jitter geom_hline labs element_blank
#' @importFrom ggplot2 scale_color_identity scale_fill_manual ggtitle
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate bind_rows
#' @importFrom reshape2 melt
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid textGrob unit gpar grid.draw
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal
#' @importFrom tidyselect all_of
#' @importFrom cowplot get_legend
#'

plot.idiffomix<-function(x, what="chromosome",CHR=1, Gene=NULL,K=3,L=3,
                         gene_cluster_name=c( "E-","E0","E+"),
                         cpg_cluster_name=c( "M-","M0","M+"),...){
    object<-x
    if(is.null(CHR)){

        warning("chromosome argument cannot be NULL for generation of plot.", call. = FALSE)
        return(invisible(NULL))
    }else{
        if(what=="chromosome"){
            tau=object$tau[[CHR]]
            pi=object$pi[[CHR]]
            temp <- matrix(tau, L, K, byrow = TRUE) * pi
            p_gene_given_cpg <- temp / rowSums(temp)
            P_gene_given_cpg <- p_gene_given_cpg
            colnames(P_gene_given_cpg) <- gene_cluster_name
            rownames(P_gene_given_cpg) <- cpg_cluster_name
            df_conditional <- melt(P_gene_given_cpg)
            colnames(df_conditional) <- c("CpG_Cluster", "Gene_Cluster", "Probability")
            df_prior <- data.frame(Cluster = as.factor( gene_cluster_name), Probability = tau)
            colnames(pi) <- gene_cluster_name
            rownames(pi) <- cpg_cluster_name
            df_pi <- melt(pi)
            colnames(df_pi) <- c("CpG_Cluster", "Gene_Cluster", "Probability")

            df_prior$Cluster <- factor(df_prior$Cluster, levels = gene_cluster_name)

            p1 <- ggplot(df_prior, aes(x=Cluster, y=Probability, fill = Probability)) +
                geom_bar(stat="identity", color="black", size=1.2)  +
                scale_fill_gradient(low="white", high="blue") +
                xlab("Gene Cluster") +
                ylim(0, 1) + theme_minimal()+ theme(legend.position = "none",
                                                    axis.text.x = element_text(face = "bold", size = 14),  # Make x-axis ticks bold
                                                    axis.text.y = element_text(face = "bold", size = 14),
                                                    axis.title.x = element_text(face = "bold", size = 14),  # Make x-axis label bold
                                                    axis.title.y = element_text(face = "bold", size = 14) )+
                geom_text(aes(label=sprintf("%.2f", Probability)),
                          vjust=-0.5,
                          size=5,              # Increase the text size
                          fontface="bold")

            df_pi$text_color <- ifelse(df_pi$Probability < 0.5, "black", "white")
            p2 <- ggplot(df_pi, aes(x=factor(Gene_Cluster), y=factor(CpG_Cluster), fill=Probability)) +
                geom_tile() +
                scale_fill_gradient(low="white", high="blue") +
                xlab("Gene Cluster") +
                ylab("P(CpG = l | Gene cluster)") +
                theme(legend.position = "none")+
                theme(
                    legend.position = "none",
                    axis.ticks.x = element_blank(), # Remove x-axis ticks
                    axis.ticks.y = element_blank(),
                    axis.text.x = element_text(face = "bold", size = 14),  # Make x-axis ticks bold
                    axis.text.y = element_text(face = "bold", size = 14),
                    axis.title.x = element_text(face = "bold", size = 14),  # Make x-axis label bold
                    axis.title.y = element_text(face = "bold", size = 14),  # Remove y-axis ticks
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.background = element_blank(), # Ensure background is clear
                    plot.background = element_blank()  )+
                geom_text(aes(label=sprintf("%.2f", Probability), color=text_color),
                          vjust=-0.5,
                          size=6,              # Increase the text size
                          fontface="bold") +
                scale_color_identity()

            # Plotting
            df_conditional$text_color <- ifelse(df_conditional$Probability < 0.5, "black", "white")
            p3 <- ggplot(df_conditional, aes(x=factor(Gene_Cluster), y=factor(CpG_Cluster), fill=Probability)) +
                geom_tile() +
                scale_fill_gradient(low="white", high="blue") +
                xlab("P(Gene = k | CpG cluster)") +
                ylab("CpG Cluster")+
                theme(
                    legend.position = "right",
                    legend.title = element_text(face = "bold"),
                    axis.ticks.x = element_blank(), # Remove x-axis ticks
                    axis.ticks.y = element_blank(), # Remove y-axis ticks
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.background = element_blank(), # Ensure background is clear
                    plot.background = element_blank(),   # Ensure the plot background is clear
                    axis.text.x = element_text(face = "bold", size = 14),  # Make x-axis ticks bold
                    axis.text.y = element_text(face = "bold", size = 14),
                    axis.title.x = element_text(face = "bold", size = 14),  # Make x-axis label bold
                    axis.title.y = element_text(face = "bold", size = 14)     # Ensure the plot background is clear
                )+
                geom_text(aes(label=sprintf("%.2f", Probability), color=text_color),
                          vjust=-0.5,
                          size=6,              # Increase the text size
                          fontface="bold") +
                scale_color_identity()

            legend <- get_legend(p3 + theme(legend.direction = "horizontal"))
            # as_ggplot(legend)

            p3 <- p3 + theme(legend.position = "none")

            # Create labels for each plot
            labelA <- textGrob("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 14))
            labelB <- textGrob("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 14))
            labelC <- textGrob("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 14))

            # Combine the plots with labels
            p1_labeled <- arrangeGrob(p1, top = labelA)
            p2_labeled <- arrangeGrob(p2, top = labelB)
            p3_labeled <- arrangeGrob(p3, top = labelC)

            # Arrange the three plots in a grid, with the common legend at the bottom
            plot_graph <- grid.arrange(
                arrangeGrob(p1_labeled, p2_labeled, p3_labeled, nrow = 1),
                legend,
                ncol = 1,
                heights = c(10, 1)  # Adjust the height ratio if needed
            )


        }else{
            final_gene_df2=object$seq_classification
            final_cpg_df2=object$meth_classification
            gene_name=Gene
            gene_col=ncol(final_gene_df2)-1
            cpg_col=ncol(final_cpg_df2)-1
            gene_reg_plot<-final_gene_df2[final_gene_df2$Gene==gene_name,]
            cpg_reg_plot<-final_cpg_df2[final_cpg_df2$Gene==gene_name,]

            max_lim=max(max(gene_reg_plot[,2:gene_col]),
                        max(cpg_reg_plot[,3:cpg_col]))+0.2
            min_lim=min(min(gene_reg_plot[,2:gene_col]),
                        min(cpg_reg_plot[,3:cpg_col]))-0.2

            columns_to_pivot<-2:gene_col
            gene_reg_plot_long <- gene_reg_plot %>%
                pivot_longer(
                    cols = all_of(columns_to_pivot),
                    names_to = "Patient",
                    values_to = "Value"
                ) %>%
                mutate(Type = "Gene Expression")


            # Prepare the CpG methylation data
            columns_to_pivot<-3:cpg_col
            cpg_reg_plot_long <- cpg_reg_plot %>%
                pivot_longer(cols = all_of(columns_to_pivot), names_to = "Patient", values_to = "Value") %>%
                mutate(Type =  CpG)

            # Combine the two datasets
            combined_df <- bind_rows(gene_reg_plot_long, cpg_reg_plot_long)
            combined_df$Final_Cluster<-ifelse(combined_df$Type=="Gene Expression","Gene",combined_df$Final_Cluster)
            combined_df$Final_Cluster <- factor(combined_df$Final_Cluster, levels = c("Gene",seq(1, L)), labels = c("Gene",cpg_cluster_name))
            gene_exp<-paste0(gene_name) #," gene")
            combined_df$Type<-ifelse(combined_df$Type=="Gene Expression",
                                     gene_exp,combined_df$Type)
            combined_df$Type <- factor(combined_df$Type, levels = unique(combined_df$Type))
            # Plotting
            p1 <- ggplot(combined_df, aes(x = Type, y = Value, fill = Final_Cluster)) +
                geom_boxplot() +
                geom_jitter(width = 0.2, size = 1.5, color = "black") +
                scale_fill_manual(values = c(
                    "Gene"="gray",
                    "M-" = "orange", "M0" = "#009E73", "M+" = "brown")) +
                xlab("") +
                ylab("Change between conditions") +
                ylim(min_lim, max_lim) +
                geom_hline(yintercept = c(0), color = "black", linewidth = 1.0, linetype = "dotted") +
                theme_minimal() +
                theme(
                    legend.position = "right",
                    axis.text.x = element_text(face = "bold", size = 16),  # Make x-axis ticks bold
                    axis.text.y = element_text(face = "bold", size = 16),
                    legend.title = element_text(face = "bold", size = 16),
                    legend.text = element_text(face = "bold", size = 16),
                    axis.title.y = element_text(face = "bold", size = 18)
                ) +
                labs(fill = "", color = "Integromix \nCluster")


            u_gk_plot<-object$U[[CHR]]
            p_gene_given_all<-unlist(u_gk_plot[u_gk_plot$Gene==gene_name,2:(K+1)])
            df_posterior <- data.frame(Cluster = as.factor(gene_cluster_name), Probability = p_gene_given_all)
            df_posterior$Cluster <- factor(df_posterior$Cluster, levels = gene_cluster_name)

            p2 <- ggplot(df_posterior, aes(x=Cluster, y=Probability,  fill = Probability)) +
                geom_bar(stat="identity", color="black", size=1.2)  +
                scale_fill_gradient(low="white", high="blue")+xlab("Gene cluster")+
                ylim(0, 1.1) +  theme_minimal()+theme(legend.position = "right",   # Ensure the plot background is clear
                                                      axis.text.x = element_text(face = "bold", size = 16),  # Make x-axis ticks bold
                                                      axis.text.y = element_text(face = "bold", size = 16),
                                                      legend.title = element_text(face = "bold", size = 14),
                                                      legend.text = element_text(face = "bold", size = 14),
                                                      axis.title.x = element_text(face = "bold", size = 16),  # Make x-axis label bold
                                                      axis.title.y = element_text(face = "bold", size = 16)   )+
                geom_text(aes(label=sprintf("%.2f", Probability)),
                          vjust=-0.5,
                          size=5,              # Increase the text size
                          fontface="bold")

            labelA <- textGrob("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))
            labelB <- textGrob("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))

            # Arrange p1 and p2 with labels
            p1_labeled <- arrangeGrob(p1, top = labelA)
            p2_labeled <- arrangeGrob(p2, top = labelB)

            # Arrange all plots into a grid
            plot_graph <- grid.arrange(
                p1_labeled,p2_labeled,
                nrow = 2
            )

        }
    }


    plot_graph

}
