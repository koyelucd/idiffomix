#' @title The idiffomix function
#' @export
#' @description Integrated differential analysis of multi-omics data using a joint mixture model
#'
#' @details This is a function to fit the joint mixture model to the transformed and filtered gene expression and methylation data.
#'
#' @param seq_data A dataframe of dimension \eqn{G \times (N+1)} containing log-fold change values for all G genes and N patients
#' @param meth_data A dataframe of dimension \eqn{C \times (N+2)} containing M-value differences between the two biological conditions for all $C$ CpG sites and $N$ patients along with the associated genes for each CpG site.
#' @param gene_chr A dataframe containing the genes and their corresponding chromosome number.
#' @param N Number of patients in the study.
#' @param K Number of clusters for expression data (default = 3; E-,E0,E+).
#' @param L Number of clusters for methylation data (default = 3; M-,M0,M+).
#' @param probs Quantile probabilities for initialization (default = c(0.1,0.9)).
#' @param parallel_process The "TRUE" option results in parallel processing of the models for increased computational efficiency. The default option has been set as "FALSE" due to package testing limitations.
#'
#' @return The function returns an object of the \code{\link[idiffomix:idiffomix]{idiffomix}} class which contains the following values:
#' \itemize{
#' \item tau - The proportion of genes belonging to K clusters.
#' \item pi - A matrix containing the probability of a CpG site belonging to cluster \emph{l}, given its associated associated gene belongs to cluster \emph{k}.
#' \item mu - The mean for each component of the gene expression data. If there is more than one component, this is a matrix whose \emph{k}th column is the mean of the \emph{k}th component of the mixture model.
#' \item sigma - The variance for each component of the gene expression data.
#' \item lambda - The mean for each component of the methylation data. If there is more than one component, this is a matrix whose \emph{l}th column is the mean of the \emph{l}th component of the mixture model.
#' \item rho - The variance for each component of the methylation data.
#' \item N - The number of patients analysed using the beta mixture models.
#' \item R - The number of sample types analysed using the beta mixture models.
#' \item U - A dataframe containing the posterior probabilities of genes belonging to the \emph{K} clusters.
#' \item V - A dataframe containing the posterior probabilities of CpG sites belonging to the \emph{L} clusters.
#' \item seq_classification - A dataframe containing the log-fold change for gene expression data and their classification corresponding to U.
#' \item meth_classification - A dataframe containing the differences in M-values for methylation data and their classification corresponding to V.
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
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom stats C quantile dnorm
#' @importFrom utils flush.console
#' @importFrom utils txtProgressBar
#' @importFrom mclust unmap
#'

idiffomix<-function(seq_data,meth_data, gene_chr,N, K=3, L=3, probs=c(0.1,0.9),
                    parallel_process = FALSE){

    seq_data_genome=seq_data
    meth_data_genome=meth_data
    chr_vec<-sort(as.numeric(unique(gene_chr$CHR)))
    u_gk_start<-list()
    v_gcl_start<-list()
    tau_start<-list()
    mu_start<-list()
    sigma_start<-list()
    pi_lk_start<-list()
    lambda_start<-list()
    rho_start<-list()
    tau_final<-list()
    mu_final<-list()
    sigma_final<-list()
    pi_lk_final<-list()
    lambda_final<-list()
    rho_final<-list()
    gene_cluster_final<-list()
    cpg_cluster_final<-list()
    u_gk_final<-list()
    v_gcl_final<-list()

    # if (nrow(seq_data_genome) == 0) {
    #     stop("seq_data_genome is empty")
    # } else {
    #     message("All checks passed: seq_data_genome and gene_chr_vec are okay.")
    # }



    initialization<-function(seq_data,meth_data,K,L,N,probs=c(0.1,0.9)){
        gene_quantile<-stats::quantile(as.vector(unlist(seq_data[,2:(N+1)])),probs = probs)
        gcluster_temp<-ifelse(rowMeans(seq_data[,2:(N+1)])<=gene_quantile[1],1,
                              ifelse(rowMeans(seq_data[,2:(N+1)])<gene_quantile[2],
                                     2,3))
        cpg_quantile<-stats::quantile(as.vector(unlist(meth_data[,3:(N+2)])),probs = probs)
        ccluster_temp<-ifelse(rowMeans(meth_data[,3:(N+2)])<=cpg_quantile[1],1,
                              ifelse(rowMeans(meth_data[,3:(N+2)])<cpg_quantile[2],
                                     2,3))


        gcluster<-gcluster_temp
        ccluster<-ccluster_temp

        ## Calculating model parameters for initialised clusters
        u_gk<-NULL
        v_gcl<-NULL

        u_gk<-mclust::unmap(gcluster, groups = 1:K)
        v_gcl<-mclust::unmap(ccluster, groups = 1:L)
        u_gk<-cbind(seq_data[,1],u_gk)
        colnames(u_gk)<-c("Gene","Cluster1","Cluster2","Cluster3")
        u_gk<-as.data.frame(u_gk)
        u_gk$Cluster1<-as.numeric(u_gk$Cluster1)
        u_gk$Cluster2<-as.numeric(u_gk$Cluster2)
        u_gk$Cluster3<-as.numeric(u_gk$Cluster3)
        v_gcl<-cbind(meth_data$Gene,meth_data$CpG,v_gcl)
        colnames(v_gcl)<-c("Gene","CpG","Cluster1","Cluster2","Cluster3")
        v_gcl<-as.data.frame(v_gcl)
        v_gcl$Cluster1<-as.numeric(v_gcl$Cluster1)
        v_gcl$Cluster2<-as.numeric(v_gcl$Cluster2)
        v_gcl$Cluster3<-as.numeric(v_gcl$Cluster3)

        return(list(u_gk=u_gk,v_gcl=v_gcl))
    }

    m_step<-function(seq_data,meth_data,u_gk,v_gcl,K,L,N,G,C){
        G=nrow(u_gk)
        C=nrow(v_gcl)
        tau<-colSums(u_gk[,2:(K+1)])/G

        pi_lk <- sapply(1:K, function(k) {
            sapply(1:L, function(l) {
                sum1 <- sum2 <- 0
                x=sapply(1:G, function(g) {
                    df_new <- v_gcl[v_gcl$Gene == u_gk$Gene[g], 3:(L+2)]
                    index <- which(v_gcl$Gene == u_gk$Gene[g])
                    Cg_num <- nrow(df_new)
                    if(Cg_num != 0)
                    {
                        sum1 <-  (u_gk[g, k + 1] * sum(df_new[, l]))
                        sum2 <-  (u_gk[g, k + 1] * Cg_num)
                    }
                    c(sum1,sum2)
                })
                sum(x[1,])/sum(x[2,])
            })
        })

        mu_vec<-vector()
        for(k in 1:K){
            temp=0
            for(g in 1:G){
                temp_seq<-as.data.frame(seq_data[seq_data$Gene==u_gk$Gene[g], -1])
                for(n in 1:N){
                    temp<-temp+(u_gk[g, k + 1] * temp_seq[,n])
                }
            }
            mu_vec[k]<-(temp)/(N*sum(u_gk[, k + 1]))
        }
        mu<-mu_vec

        sigma_sq<-vector()
        for(k in 1:K){
            temp=0
            for(g in 1:G){
                temp_seq<-as.data.frame(seq_data[seq_data$Gene==u_gk$Gene[g], -1])
                for(n in 1:N){
                    temp<-temp+(u_gk[g, k + 1] * ((temp_seq[,n] - mu_vec[k])^2))
                }
            }
            sigma_sq[k]<-(temp)/(N*sum(u_gk[, k + 1]))
        }
        sigma<-sigma_sq^0.5

        lambda_vec <- sapply(1:L, function(l) {
            sum1 <- sum(sapply(1:C, function(c) {
                sum(v_gcl[c, l + 2] * meth_data[meth_data$CpG==v_gcl$CpG[c], 3:(N+2)])
            }))
            sum1 / (N * sum(v_gcl[, l + 2]))
        })
        lambda<-lambda_vec

        rho_sq <- numeric(L)
        for (l in 1:L) {
            Gc<-unique(v_gcl$Gene)
            sum2 <- numeric(length(Gc))
            sum4<-numeric(length(Gc))
            for (g in 1:length(Gc)) {
                cpg <- meth_data[meth_data$Gene == Gc[g], ]
                v_gcl_temp <- v_gcl[v_gcl$Gene == Gc[g], ]
                sum1 <- numeric(nrow(cpg))
                for (c in 1:nrow(cpg)) {
                    temp <- numeric(N)
                    for (n in 1:N) {
                        temp[n] <- v_gcl_temp[c, l + 2] *
                            ((cpg[cpg$CpG==v_gcl_temp$CpG[c], n + 2] -
                                  lambda_vec[l])^2)
                    }
                    sum1[c] <- sum(temp)
                }
                sum4[g] <- sum(v_gcl_temp[, l + 2])
                sum2[g] <- sum(sum1)
            }
            rho_sq[l] <- sum(sum2) / (N * sum(sum4))
        }
        rho<-rho_sq^0.5


        ## Constraining the standard deviations to be same for each cluster
        sigma_temp<-sum(tau*sigma)
        sigma<-rep(sigma_temp,times=K)
        pi_temp<-rowSums(t(pi_lk)*(tau))
        rho_temp<-sum(pi_temp*rho)
        rho<-rep(rho_temp,times=L)

        return(list(tau=tau,pi_lk=pi_lk,mu=mu,sigma=sigma,
                    lambda=lambda,rho=rho))
    }
    e_step<-function(seq_data,meth_data,u_gk,v_gcl,K,L,N,G,C,tau,pi_lk,mu,sigma,lambda,rho,conv_th=0.001){

        ugk_vgcl_s<-list()
        u_gk_s<-list()
        v_gcl_s<-list()

        has_converged=FALSE
        convergence_threshold=conv_th
        change_measure<-vector()
        change_measure2<-vector()
        s=1
        while(!has_converged){
            u_gk_old<-u_gk
            v_gcl_old<-v_gcl

            ####### u_gk estimation ##############################
            prob_u_gk <- matrix(0, nrow = G, ncol = K+1)
            for(g in 1:G){
                for(k in 1:K){
                    temp1<-sum(sapply(1:N, function(n)
                        (stats::dnorm(seq_data[g, n + 1], mu[k],
                               sigma[k],log=T))))
                    index<-which(v_gcl$Gene==seq_data$Gene[g])
                    Cg_num=length(index)
                    temp2<-0
                    if(Cg_num!=0)
                    {
                        temp2 <-
                            sum(sapply(1:Cg_num, function(c) {
                                sum(sapply(1:L, function(l) {
                                    if(v_gcl[index[c], l+2]!=0 & pi_lk[l,k]!=0)
                                    {
                                        v_gcl[index[c], l+2]*
                                            log(pi_lk[l, k])
                                    }else {0}
                                }))
                            }))
                    }
                    prob_u_gk[g,1]<-seq_data$Gene[g]
                    prob_u_gk[g,k+1]<-
                        log(tau[k])+
                        temp1+
                        temp2
                }
            }
            colnames(prob_u_gk)<-colnames(u_gk)
            prob_u_gk<-as.data.frame(prob_u_gk)
            cols_to_convert <- c("Cluster1", "Cluster2","Cluster3")
            prob_u_gk[cols_to_convert] <- lapply(prob_u_gk[cols_to_convert], as.numeric)
            temp3=apply(prob_u_gk[,2:(K+1)],1,max)
            prob_u_gk[,2:(K+1)]=exp(prob_u_gk[,2:(K+1)]-(temp3))
            prob_u_gk[,2:(K+1)]<-prob_u_gk[,2:(K+1)]/rowSums(prob_u_gk[,2:(K+1)])
            u_gk_s[[s]]<-prob_u_gk
            g_final2<-apply(prob_u_gk[,2:(K+1)],1,which.max)
            u_gk[,2:(K+1)]<-mclust::unmap(g_final2,group=1:K)


            # ####### v_gcl estimation ##############################
            prob_v_gcl<-matrix(0,nrow=C,ncol=L+2)
            for(c in 1:C){
                for(l in 1:L){
                    temp1 <- sum(sapply(1:N, function(n)
                        (stats::dnorm(meth_data[c,n+2],lambda[l],
                               rho[l],log = T))))
                    index<-which(u_gk$Gene==meth_data$Gene[c])
                    temp2<-0
                    if(length(index)!=0)
                    {
                        temp2 <- sum(sapply(1:K, function(k) {
                            if(u_gk[index, k+1]!=0 & pi_lk[l,k]!=0)
                            {
                                u_gk[index, k+1]*log(pi_lk[l, k])
                            }else {0}
                        }))
                    }

                    prob_v_gcl[c,l+2]<-temp1+temp2
                }
                prob_v_gcl[c,1]<-meth_data[c,1]
                prob_v_gcl[c,2]<-meth_data[c,2]
            }
            colnames(prob_v_gcl)<-colnames(v_gcl)
            prob_v_gcl<-as.data.frame(prob_v_gcl)
            cols_to_convert <- c("Cluster1", "Cluster2","Cluster3")
            prob_v_gcl[cols_to_convert] <- lapply(prob_v_gcl[cols_to_convert], as.numeric)
            temp3=apply(prob_v_gcl[,3:(L+2)],1,max)
            prob_v_gcl[,3:(L+2)]=exp(prob_v_gcl[,3:(L+2)]-(temp3))
            prob_v_gcl[,3:(L+2)]<-prob_v_gcl[,3:(L+2)]/rowSums(prob_v_gcl[,3:(L+2)])
            v_gcl_s[[s]]<-prob_v_gcl
            c_final2<-apply(prob_v_gcl[,3:(L+2)],1,which.max)
            v_gcl[,3:(L+2)]<-mclust::unmap(c_final2,group=1:L)


            if(s>1){
                change_measure[s] <- (mean(as.matrix(abs(u_gk[,2:(K+1)] - u_gk_old[,2:(K+1)]),na.rm=TRUE)))
                change_measure2[s] <- (mean(as.matrix(abs(v_gcl[,3:(L+2)] - v_gcl_old[,3:(L+2)]))))
                if (change_measure[s] < convergence_threshold &
                    change_measure2[s] < convergence_threshold) {
                    has_converged <- TRUE
                }
            }else{
                change_measure[s]<-0.002
                change_measure2[s]<-0.002
                has_converged<-FALSE
            }

            s=s+1
        }
        len_s=length(u_gk_s)
        u_gk<-u_gk_s[[len_s]]
        u_gk<-u_gk[order(u_gk$Gene),]
        u_gk[,2:(K+1)]<-u_gk[,2:(K+1)]/length(u_gk_s)
        u_gk[,2:(K+1)]<-u_gk[,2:(K+1)]/rowSums(u_gk[,2:(K+1)])

        v_gcl<-v_gcl_s[[len_s]]
        v_gcl<-v_gcl[order(v_gcl$Gene,v_gcl$CpG),]
        v_gcl[,3:(L+2)]<-v_gcl[,3:(L+2)]/length(v_gcl_s)
        v_gcl[,3:(L+2)]<-v_gcl[,3:(L+2)]/rowSums(v_gcl[,3:(L+2)])

        return(list(u_gk=u_gk,v_gcl=v_gcl))
    }

    if (parallel_process == FALSE) {
        # use 2 cores in CRAN/Travis/AppVeyor
        ncores <- 2L
    } else {
        # use all cores in devtools::test()
        ncores <- parallel::detectCores()
    }
    num_cores <- ncores - 1  # Adjust as needed
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    tau_final_par<-list()
    mu_final_par<-list()
    sigma_final_par<-list()
    pi_lk_final_par<-list()
    lambda_final_par<-list()
    rho_final_par<-list()
    gene_cluster_final_par<-list()
    cpg_cluster_final_par<-list()
    u_gk_final_par<-list()
    v_gcl_final_par<-list()

    results <- foreach(chr = chr_vec, .combine = list, .multicombine = TRUE, .packages = c('mclust')) %dopar% {

        ## Select data for the corresponding CHR
        gene_chr_vec <- as.vector(unlist(gene_chr[gene_chr$CHR == chr, "Gene"]))
        seq_transformed_df <- seq_data_genome[seq_data_genome$Gene %in% gene_chr_vec, ]
        meth_transformed_df <- meth_data_genome[meth_data_genome$Gene %in% gene_chr_vec, ]
        seq_transformed_df <- seq_transformed_df[order(seq_transformed_df$Gene), ]
        meth_transformed_df <- meth_transformed_df[order(meth_transformed_df$Gene, meth_transformed_df$CpG), ]

        ####################### Coordinate ascent method #########################
        ############## Initialization #######################
        G = nrow(seq_transformed_df)
        C = nrow(meth_transformed_df)
        iter = 20
        epsilon = 1e-5

        ## Quantile clustering
        initial_clustering <- initialization(seq_transformed_df, meth_transformed_df, K, L, N, probs = probs)
        u_gk <- initial_clustering$u_gk
        v_gcl <- initial_clustering$v_gcl

        ########## M-step for estimating the parameters ################
        initial_output <- m_step(seq_transformed_df, meth_transformed_df, u_gk, v_gcl, K, L, N,G,C)

        ## Store the initialized parameters
        tau <- initial_output$tau
        mu <- initial_output$mu
        sigma <- initial_output$sigma
        pi_lk <- initial_output$pi_lk
        lambda <- initial_output$lambda
        rho <- initial_output$rho

        total_diff <- vector()

        ################# EM algorithm ########################
        for (i in 1:100) {
            ######### E-Step #####################
            e_step_output <- e_step(seq_transformed_df, meth_transformed_df, u_gk, v_gcl, K, L, N,G,C, tau, pi_lk, mu, sigma, lambda, rho)

            ## Update u_gk and v_gcl
            u_gk <- e_step_output$u_gk
            v_gcl <- e_step_output$v_gcl

            u_gk_prob <- e_step_output$u_gk
            v_gcl_prob <- e_step_output$v_gcl

            ################## M-step #########################
            g_final2 <- apply(u_gk[, 2:(K+1)], 1, which.max)
            u_gk[,2:(K+1)] <- as.data.frame(mclust::unmap(g_final2, group = 1:K))
            colnames(u_gk) <- c("Gene", "Cluster1", "Cluster2","Cluster3")

            c_final2 <- apply(v_gcl[, 3:(L+2)], 1, which.max)
            v_gcl[,3:(L+2)] <- as.data.frame(mclust::unmap(c_final2, group = 1:L))
            colnames(v_gcl) <- c("Gene", "CpG", "Cluster1", "Cluster2","Cluster3")

            m_step_output <- m_step(seq_transformed_df, meth_transformed_df, u_gk, v_gcl, K, L, N)

            ## Store the new estimates
            tau_est <- m_step_output$tau
            pi_lk_est <- m_step_output$pi_lk
            mu_est <- m_step_output$mu
            sigma_est <- m_step_output$sigma
            lambda_est <- m_step_output$lambda
            rho_est <- m_step_output$rho

            ## Calculate convergence criteria
            pi_diff <- mean(pi_lk_est - pi_lk)
            tau_diff <- mean(tau_est - tau)
            mu_diff <- mean(mu_est - mu)
            sigma_diff <- mean(sigma_est - sigma)
            lambda_diff <- mean(lambda_est - lambda)
            rho_diff <- mean(rho_est - rho)

            total_diff[i] <- (pi_diff + tau_diff + mu_diff + sigma_diff + lambda_diff + rho_diff) / 6

            ## Update parameters
            pi_lk <- pi_lk_est
            tau <- tau_est
            mu <- mu_est
            sigma <- sigma_est
            lambda <- lambda_est
            rho <- rho_est

            if (abs(total_diff[i]) < epsilon) {
                break
            }
        }
        gene_cluster_final_est = as.factor(apply(u_gk_prob[, 2:(K+1)], 1, which.max))
        cpg_cluster_final_est = as.factor(apply(v_gcl_prob[, 3:(L+2)], 1, which.max))
        ## Return results for the current chromosome
        list(
            tau_final_est = tau_est,
            mu_final_est = mu_est,
            sigma_final_est = sigma_est,
            pi_lk_final_est = pi_lk_est,
            lambda_final_est = lambda_est,
            rho_final_est = rho_est,
            gene_cluster_final_est = gene_cluster_final_est,
            cpg_cluster_final_est = cpg_cluster_final_est,
            u_gk_final_est = u_gk_prob,
            v_gcl_final_est = v_gcl_prob
        )
    }

    # Stop the cluster after use
    stopCluster(cl)


    for (i in seq_along(chr_vec)) {
        chr <- chr_vec[i]

        # Access the individual results for each chromosome
        tau_final_par[[chr]] <- results[[i]]$tau_final_est
        mu_final_par[[chr]] <- results[[i]]$mu_final_est
        sigma_final_par[[chr]] <- results[[i]]$sigma_final_est
        pi_lk_final_par[[chr]] <- results[[i]]$pi_lk_final_est
        lambda_final_par[[chr]] <- results[[i]]$lambda_final_est
        rho_final_par[[chr]] <- results[[i]]$rho_final_est
        gene_cluster_final_par[[chr]] <- results[[i]]$gene_cluster_final_est
        cpg_cluster_final_par[[chr]] <- results[[i]]$cpg_cluster_final_est
        u_gk_final_par[[chr]] <- results[[i]]$u_gk_final_est
        v_gcl_final_par[[chr]] <- results[[i]]$v_gcl_final_est
    }

    final_gene_df_par<-as.data.frame(matrix(NA,nrow=0,ncol=ncol(seq_data_genome)+1))
    colnames(final_gene_df_par)<-c(colnames(seq_data_genome),"Final_Cluster")
    final_cpg_df_par<-as.data.frame(matrix(NA,nrow=0,ncol=ncol(meth_data_genome)+1))
    colnames(final_cpg_df_par)<-c(colnames(meth_data_genome),"Final_Cluster")
    for(chr in chr_vec){
        # print(chr)

        gene_chr_vec<-as.vector(unlist(gene_chr[gene_chr$CHR==chr,"Gene"]))
        seq_transformed_df_2<-seq_data_genome[seq_data_genome$Gene %in% gene_chr_vec,]
        meth_transformed_df2<-meth_data_genome[meth_data_genome$Gene %in% gene_chr_vec,]
        seq_transformed_df_2<-seq_transformed_df_2[order(seq_transformed_df_2$Gene),]
        meth_transformed_df2<-meth_transformed_df2[order(meth_transformed_df2$Gene,
                                                         meth_transformed_df2$CpG),]
        seq_transformed_df_2$Final_Cluster<-gene_cluster_final_par[[chr]]
        meth_transformed_df2$Final_Cluster<-cpg_cluster_final_par[[chr]]
        final_gene_df_par<-rbind(final_gene_df_par,seq_transformed_df_2)
        final_cpg_df_par<-rbind(final_cpg_df_par,meth_transformed_df2)

    }

    idiffomix_out<-list(tau=tau_final_par,
                        pi=pi_lk_final_par,
                        mu=mu_final_par,
                        sigma=sigma_final_par,
                        lambda=lambda_final_par,
                        rho=rho_final_par,
                        U=u_gk_final_par,
                        V=v_gcl_final_par,
                        seq_classification = final_gene_df_par,
                        meth_classification = final_cpg_df_par)
    class(idiffomix_out)<-"idiffomix"
    return(idiffomix_out)

}
