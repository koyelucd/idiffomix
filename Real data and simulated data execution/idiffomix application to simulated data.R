install.packages("fourPNO")
library(fourPNO)    
library(mclust)
library(readr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(limma)
library(purrr)
library(edgeR)
library(e1071)
library(MASS)
### Simulate data
perf_metrics<-matrix(NA,nrow=100,ncol=4)
perf_metrics_gene_intmod<-as.data.frame(perf_metrics)
colnames(perf_metrics_gene_intmod)<-c("FDR","Sens","Spec","ARI")
perf_metrics_gene_mclust<-as.data.frame(perf_metrics)
colnames(perf_metrics_gene_mclust)<-c("FDR","Sens","Spec","ARI")
perf_metrics_cpg_intmod<-as.data.frame(perf_metrics)
colnames(perf_metrics_cpg_intmod)<-c("FDR","Sens","Spec","ARI")
perf_metrics_cpg_mclust<-as.data.frame(perf_metrics)
colnames(perf_metrics_cpg_mclust)<-c("FDR","Sens","Spec","ARI")
perf_metrics_gene_limma<-as.data.frame(perf_metrics)
colnames(perf_metrics_gene_limma)<-c("FDR","Sens","Spec","ARI")
perf_metrics_cpg_limma<-as.data.frame(perf_metrics)
colnames(perf_metrics_cpg_limma)<-c("FDR","Sens","Spec","ARI")

### case 3 initialization
# perf_metrics<-matrix(NA,nrow=100,ncol=2)
# perf_metrics_gene_intmod<-as.data.frame(perf_metrics)
# colnames(perf_metrics_gene_intmod)<-c("ARI_mclust","ARI_limma")
# perf_metrics_cpg_intmod<-as.data.frame(perf_metrics)
# colnames(perf_metrics_cpg_intmod)<-c("ARI_mclust","ARI_limma")
# perf_metrics_cpg_mclust<-as.data.frame(perf_metrics)



seed_val<-round(runif(100,min=0,max=1)*1000)


df_gene_list<-list()
df_cpg_list<-list()
Cg_total<-0

for(seed in 3:100){
    print(seed)
    seed=56
    set.seed(seed_val[seed])
    K_0<-3
    G<-500
    N<-4
    tau_0<-c(0.1,0.8,0.1)
    mu_0<-c(-1.5,0,1.5)
    sigma_0<-1
    index_g<-sample(1:K_0,size=G,prob=tau_0,replace = TRUE)
    index_non<-sample(1:2,size=table(index_g)[2],replace = TRUE)
    df_gene<-matrix(NA,nrow=G,ncol=(N*2))
    generate_row <- function(ind) {
        if(ind==3){
            c(rnegbin(4,mu=10000,theta = 5),rnegbin(4,mu=60000,theta = 5))
        }else if(ind ==1){
                c(rnegbin(4,mu=10000,theta = 5),rnegbin(4,mu=4000,theta = 5))
        }else{
            c(rnegbin(4,mu=10000,theta = 5),rnegbin(4,mu=10000,theta = 5))
            
        }
    }
    
    
    df_gene <- t(sapply(index_g, generate_row))
    gene_name<-paste("gene_",seq(1:G),sep = "")
    df_gene<-cbind(gene_name,df_gene,index_g)
    df_gene<-as.data.frame(df_gene)
    colnames(df_gene)<-c("Gene_name","Patient1_GX1","Patient2_GX1","Patient3_GX1","Patient4_GX1","Patient1_GX2","Patient2_GX2","Patient3_GX2","Patient4_GX2","Group_Gene")
    df_gene$Patient1_GX1<-as.numeric(df_gene$Patient1_GX1)
    df_gene$Patient2_GX1<-as.numeric(df_gene$Patient2_GX1)
    df_gene$Patient3_GX1<-as.numeric(df_gene$Patient3_GX1)
    df_gene$Patient4_GX1<-as.numeric(df_gene$Patient4_GX1)
    df_gene$Patient1_GX2<-as.numeric(df_gene$Patient1_GX2)
    df_gene$Patient2_GX2<-as.numeric(df_gene$Patient2_GX2)
    df_gene$Patient3_GX2<-as.numeric(df_gene$Patient3_GX2)
    df_gene$Patient4_GX2<-as.numeric(df_gene$Patient4_GX2)
    df_gene$Group_Gene<-as.factor(df_gene$Group_Gene)
    add_poisson_noise <- function(x) {
        rpois(1, lambda = x)
    }

    # Apply the function to each element of the data frame
    noisy_data <- as.data.frame(lapply(df_gene[,c(2:9)], function(column) {
        sapply(column, function(x) add_poisson_noise(x))
    }))
    df_gene[,2:9]<-noisy_data
    
    

    
    # View the first few rows of the noisy data
    df_gene_raw<-df_gene
    df_gene_list[[seed]]<-df_gene_raw
    rownames(df_gene)<-df_gene$Gene_name
    gene_df_new<-df_gene[,-1]
    gene_df_new<-gene_df_new[,-9]
    
    targets <- data.frame(Patients=c(1,2,3,4,1,2,3,4),
                          Conditions=c("Benign","Benign",
                                       "Benign","Benign","Tumour",
                                       "Tumour","Tumour","Tumour"))
    Patients <- factor(targets$Patients)
    Conditions <- factor(targets$Conditions, levels=c("Benign","Tumour"))
    design <- model.matrix(~0+Patients+Conditions)
    design
    
    keep <- filterByExpr(gene_df_new, design)
    gene_df_new2 <- gene_df_new[keep,]
    gene_cpm<-gene_df_new2
    # Create DGEList object
    dge <- DGEList(counts = gene_cpm)
    nonzero <- rowSums(dge$counts) > 0
    dge %>% .[nonzero,]
    dge %>% calcNormFactors
    cpm_values <- cpm(dge)
    
    log_transformed<-log2(cpm_values + 1)
    # log_transformed<-log2(gene_df_new+1)
    log_transformed<-as.data.frame(log_transformed)
    log_transformed<-as.data.frame(cbind(rownames(log_transformed),log_transformed))
    log_transformed[,2]<-log_transformed[,6]-log_transformed[,2]
    log_transformed[,3]<-log_transformed[,7]-log_transformed[,3]
    log_transformed[,4]<-log_transformed[,8]-log_transformed[,4]
    log_transformed[,5]<-log_transformed[,9]-log_transformed[,5]
    log_transformed<-log_transformed[,-(6:9)]
    df_gene<-log_transformed
    df_gene<-as.data.frame(cbind(df_gene,df_gene_raw$Group_Gene))
    colnames(df_gene)<-c("Gene","Patient1_GX","Patient2_GX","Patient3_GX","Patient4_GX","Group_Gene")
    
    
    df_gene<-df_gene[order(df_gene$Gene),]
    L_0 <- 3
    #### Different pi_lk settings
    # pi_0<-matrix(c(0.1,0.05,0.4,0.5,0.90,0.5,0.4,0.05,0.1), L_0, L_0,byrow = T) ##case 1
    # pi_0<-matrix(c(0.1,0.1,0.8,0.1,0.80,0.1,0.8,0.1,0.1), L_0, L_0,byrow = T)   ##case 2
    pi_0<-matrix(c(0.2,0.2,0.2,0.6,0.6,0.6,0.2,0.2,0.2), L_0, L_0,byrow = T)      ##case 3
    lambda_0<-c(-1.5,0,1.5)
    rho_0<-1
    Cg<-floor(runif(500,min=0.1,max=1)*30)
    Cg_total<-Cg_total+Cg
    cpg_ind <- rep(NA, sum(Cg))
    cumsum_Cg<-cumsum(Cg)
    start=1
    end=cumsum_Cg[1]
    for (g in 1:G)
    {
        # print(g)
        cpg_ind[start:end] <- sample(1:L_0, size = Cg[g], prob = pi_0[, df_gene$Group_Gene[g]], replace = TRUE)
        if(g!=G){
            start=end+1
            end=cumsum_Cg[g+1]
        }
        
    }
    generate_row <- function(ind) {
        rmvnorm(1, rep(lambda_0[ind], N), diag(rho_0, N))
        if(ind==1){
            c(rbeta(4,20,3),rbeta(4,3,20))
        }else if(ind ==3){
            c(rbeta(4,3,20),rbeta(4,20,3))
        }else{
            index_non=sample(1:3,size=1,replace = TRUE)
            if(index_non==3){
                c(rbeta(4,20,3),rbeta(4,20,3))
            }else if(index_non==1){
                c(rbeta(4,3,20),rbeta(4,3,20))
            }else{
                c(rbeta(4,4,3),rbeta(4,4,3))
            }
        }
    }
    df_cpg <- t(sapply(cpg_ind, generate_row))
    cpg_name=NULL
    for(g in 1:G){
        for(c in 1:Cg[g]){
            x<-unlist(strsplit(df_gene$Gene[g],"_"))[2]
            cpg_name<-c(cpg_name,paste("CpG",x,c,sep="_"))
        }
    }
    
    df_cpg<-cbind(cpg_name,df_cpg,cpg_ind)
    df_cpg<-as.data.frame(df_cpg)
    colnames(df_cpg)<-c("CpG","Patient1_M1","Patient2_M1","Patient3_M1","Patient4_M1","Patient1_M2","Patient2_M2","Patient3_M2","Patient4_M2","Group_CpG")
    df_cpg$Patient1_M1<-as.numeric(df_cpg$Patient1_M1)
    df_cpg$Patient2_M1<-as.numeric(df_cpg$Patient2_M1)
    df_cpg$Patient3_M1<-as.numeric(df_cpg$Patient3_M1)
    df_cpg$Patient4_M1<-as.numeric(df_cpg$Patient4_M1)
    df_cpg$Patient1_M2<-as.numeric(df_cpg$Patient1_M2)
    df_cpg$Patient2_M2<-as.numeric(df_cpg$Patient2_M2)
    df_cpg$Patient3_M2<-as.numeric(df_cpg$Patient3_M2)
    df_cpg$Patient4_M2<-as.numeric(df_cpg$Patient4_M2)
    df_cpg$Group_CpG<-as.factor(df_cpg$Group_CpG)
    min_x<-min(df_cpg[,2:9])
    max_x<-max(df_cpg[,2:9])
    if(min_x<=0){
        values <- as.vector(as.matrix(df_cpg[, 2:9]))
        unique_values <- sort(unique(values))
        min_x <- unique_values[2]
    }
    if(max_x>=1){
        values <- as.vector(as.matrix(df_cpg[, 2:9]))
        unique_values <- sort(unique(values))
        max_x <- unique_values[length(unique_values) - 1]
    }
    add_gaussian_noise <- function(x) {
       x= x + rnorm(1, mean = 0, sd = 0.05)
        if(x<=0){
            x=min_x
        }else if(x>=1){
            x=max_x
        }
       return(x)
    }
    
    # Apply the function to each element of the dataframe
    noisy_beta_data <- as.data.frame(lapply(df_cpg[,2:9], function(column) {
        sapply(column, add_gaussian_noise)
    }))
    df_cpg[,2:9]<-noisy_beta_data
    
    new_df <- as.data.frame(df_gene[rep(1:nrow(df_gene), times = Cg), ])
    new_df<-new_df[order(new_df$Gene),]
    df_complete<-cbind(new_df$Gene,df_cpg,new_df$Group_Gene)
    df_complete_raw<-df_complete
    df_cpg_list[[seed]]<-df_complete_raw
    colnames(df_complete)[1]<-c("Gene")
    colnames(df_complete)[12]<-"Group_Gene"
    log_transformed<-df_complete
    log_transformed[,3:10]<-log(log_transformed[,3:10]/(1-log_transformed[,3:10]))
    log_transformed[,3]<-log_transformed[,7]-log_transformed[,3]
    log_transformed[,4]<-log_transformed[,8]-log_transformed[,4]
    log_transformed[,5]<-log_transformed[,9]-log_transformed[,5]
    log_transformed[,6]<-log_transformed[,10]-log_transformed[,6]
    log_transformed<-log_transformed[,-(7:10)]
    df_complete<-log_transformed
    colnames(df_complete)<-c("Gene","CpG","Patient1_M","Patient2_M","Patient3_M","Patient4_M","Group_CpG","Group_Gene")
    
    df_cpg<-df_complete[,1:6]
    plot(density(df_cpg$Patient1_M))
    
    pi_l1<-table(df_complete[df_complete$Group_Gene==1,7])/nrow(df_complete[df_complete$Group_Gene==1,])
    pi_l2<-table(df_complete[df_complete$Group_Gene==2,7])/nrow(df_complete[df_complete$Group_Gene==2,])
    pi_l3<-table(df_complete[df_complete$Group_Gene==3,7])/nrow(df_complete[df_complete$Group_Gene==3,])
    pi_lk_0<-cbind(pi_l1,pi_l2,pi_l3)
    pi_lk_0
    
    
    seq_transformed_df<-df_gene[,1:5]
    meth_transformed_df<-df_complete[,c(1:6)]
    K=3
    L=3
    C=nrow(meth_transformed_df)

    
    ####################### Coordinate ascent method #########################
    ############## Initialization #######################
    ## Quantile clustering
    gene_quantile<-quantile(as.vector(unlist(seq_transformed_df[,2:5])),probs = c(0.1,0.90))
    updated_gcluster<-ifelse(rowMeans(seq_transformed_df[,2:5])<=gene_quantile[1],1,ifelse(rowMeans(seq_transformed_df[,2:5])<gene_quantile[2],2,3))
    cpg_quantile<-quantile(as.vector(unlist(meth_transformed_df[,3:6])),probs = c(0.15,0.85))
    updated_ccluster<-ifelse(rowMeans(meth_transformed_df[,3:6])<=cpg_quantile[1],1,ifelse(rowMeans(meth_transformed_df[,3:6])<cpg_quantile[2],2,3))
    
    
    
    u_gk_start<-unmap(updated_gcluster, groups = 1:K)
    v_gcl_start<-unmap(updated_ccluster, groups = 1:L)
    u_gk_start<-cbind(seq_transformed_df$Gene,u_gk_start)
    colnames(u_gk_start)<-c("Gene","Cluster1","Cluster2","Cluster3")
    u_gk_start<-as.data.frame(u_gk_start)
    u_gk_start$Cluster1<-as.numeric(u_gk_start$Cluster1)
    u_gk_start$Cluster2<-as.numeric(u_gk_start$Cluster2)
    u_gk_start$Cluster3<-as.numeric(u_gk_start$Cluster3)
    v_gcl_start<-cbind(meth_transformed_df$Gene,meth_transformed_df$CpG,v_gcl_start)
    colnames(v_gcl_start)<-c("Gene","CpG","Cluster1","Cluster2","Cluster3")
    v_gcl_start<-as.data.frame(v_gcl_start)
    v_gcl_start$Cluster1<-as.numeric(v_gcl_start$Cluster1)
    v_gcl_start$Cluster2<-as.numeric(v_gcl_start$Cluster2)
    v_gcl_start$Cluster3<-as.numeric(v_gcl_start$Cluster3)

    
    ########## M-step for estimating the parameters ################
    tau_start<-colSums(u_gk_start[,2:(K+1)])/G
    
    pi_lk_start <- sapply(1:K, function(k) {
        sapply(1:L, function(l) {
            sum1 <- sum2 <- 0
            x=sapply(1:G, function(g) {
                df_new <- v_gcl_start[v_gcl_start$Gene == u_gk_start$Gene[g], 3:(L+2)]
                index <- which(v_gcl_start$Gene == u_gk_start$Gene[g])
                Cg_num <- nrow(df_new)
                if(Cg_num != 0)
                {
                    sum1 <-  (u_gk_start[g, k + 1] * sum(df_new[, l]))
                    sum2 <-  (u_gk_start[g, k + 1] * Cg_num)  
                }
                c(sum1,sum2)
            })
            sum(x[1,])/sum(x[2,])
        })
    })
    
    mu_start_vec<-vector()
    for(k in 1:K){
        temp=0
        for(g in 1:G){
            temp_seq<-as.data.frame(seq_transformed_df[seq_transformed_df$Gene==u_gk_start$Gene[g], -1])
            for(n in 1:N){
                temp<-temp+(u_gk_start[g, k + 1] * temp_seq[,n])
            }
        }
        mu_start_vec[k]<-(temp)/(N*sum(u_gk_start[, k + 1]))
    }
    mu_start<-mu_start_vec
    
    sigma_sq_start<-vector()
    for(k in 1:K){
        temp=0
        for(g in 1:G){
            temp_seq<-as.data.frame(seq_transformed_df[seq_transformed_df$Gene==u_gk_start$Gene[g], -1])
            for(n in 1:N){
                temp<-temp+(u_gk_start[g, k + 1] * ((temp_seq[,n] - mu_start_vec[k])^2))
            }
        }
        sigma_sq_start[k]<-(temp)/(N*sum(u_gk_start[, k + 1]))
    }
    sigma_start<-sigma_sq_start^0.5
    
    lambda_start_vec <- sapply(1:L, function(l) {
        sum1 <- sum(sapply(1:C, function(c) {
            sum(v_gcl_start[c, l + 2] * meth_transformed_df[meth_transformed_df$CpG==v_gcl_start$CpG[c], 3:(N+2)])
        }))
        sum1 / (N * sum(v_gcl_start[, l + 2]))
    })
    lambda_start<-lambda_start_vec
    
    rho_sq_start <- numeric(L)
    for (l in 1:L) {
        Gc<-unique(v_gcl_start$Gene)
        sum2 <- numeric(length(Gc))
        sum4<-numeric(length(Gc))
        for (g in 1:length(Gc)) {
            cpg <- meth_transformed_df[meth_transformed_df$Gene == Gc[g], ]
            v_gcl_temp <- v_gcl_start[v_gcl_start$Gene == Gc[g], ]
            sum1 <- numeric(nrow(cpg))
            for (c in 1:nrow(cpg)) {
                temp <- numeric(N)
                for (n in 1:N) {
                    temp[n] <- v_gcl_temp[c, l + 2] * ((cpg[cpg$CpG==v_gcl_temp$CpG[c], n + 2] - lambda_start_vec[l])^2)
                }
                sum1[c] <- sum(temp)
            }
            sum4[g] <- sum(v_gcl_temp[, l + 2])
            sum2[g] <- sum(sum1)
        }
        rho_sq_start[l] <- sum(sum2) / (N * sum(sum4))
    }
    rho_start<-rho_sq_start^0.5
    
    
    ## Constraining the standard deviations to be same for each cluster
    sigma_temp<-sum(tau_start*sigma_start)
    sigma_start<-rep(sigma_temp,times=K)
    pi_temp<-rowSums(t(pi_lk_start)*(tau_start))
    rho_temp<-sum(pi_temp*rho_start)
    rho_start<-rep(rho_temp,times=L)
    
    u_gk<-u_gk_start #[,2:(K+1)]
    v_gcl<-v_gcl_start #[,3:(L+2)]
    tau<-tau_start
    mu<-mu_start
    sigma<-sigma_start
    pi_lk<-pi_lk_start
    lambda<-lambda_start
    rho<-rho_start
    total_diff=vector()
    u_gk_exp<-u_gk
    v_gcl_exp<-v_gcl

    epsilon=1e-5
    iter=10
    ################# EM algorithm ########################
    for(i in 1:100){
        ################ E-step ###################################
        ugk_vgcl_s<-list()    
        u_gk_s<-list()
        v_gcl_s<-list()
        
        has_converged=FALSE
        convergence_threshold=0.001
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
                        (dnorm(seq_transformed_df[g, n + 1], mu[k],
                               sigma[k],log=T))))
                    index<-which(v_gcl$Gene==seq_transformed_df$Gene[g])
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
                    prob_u_gk[g,1]<-seq_transformed_df$Gene[g]
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
            temp3=apply(prob_u_gk[,2:4],1,max)
            prob_u_gk[,2:4]=exp(prob_u_gk[,2:4]-(temp3))
            prob_u_gk[,2:4]<-prob_u_gk[,2:4]/rowSums(prob_u_gk[,2:4])
            u_gk_s[[s]]<-prob_u_gk
            g_final2<-apply(prob_u_gk[,2:4],1,which.max)
            u_gk[,2:4]<-unmap(g_final2,group=1:K)
            
            
            # ####### v_gcl estimation ##############################
            prob_v_gcl<-matrix(0,nrow=C,ncol=L+2)
            for(c in 1:C){
                for(l in 1:L){
                    temp1 <- sum(sapply(1:N, function(n)
                        (dnorm(meth_transformed_df[c,n+2],lambda[l], 
                               rho[l],log = T))))
                    index<-which(u_gk$Gene==meth_transformed_df$Gene[c])
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
                prob_v_gcl[c,1]<-meth_transformed_df[c,1]
                prob_v_gcl[c,2]<-meth_transformed_df[c,2]
            }
            colnames(prob_v_gcl)<-colnames(v_gcl)
            prob_v_gcl<-as.data.frame(prob_v_gcl)
            cols_to_convert <- c("Cluster1", "Cluster2","Cluster3")
            prob_v_gcl[cols_to_convert] <- lapply(prob_v_gcl[cols_to_convert], as.numeric)
            temp3=apply(prob_v_gcl[,3:5],1,max)
            prob_v_gcl[,3:5]=exp(prob_v_gcl[,3:5]-(temp3))
            prob_v_gcl[,3:5]<-prob_v_gcl[,3:5]/rowSums(prob_v_gcl[,3:5])
            v_gcl_s[[s]]<-prob_v_gcl
            c_final2<-apply(prob_v_gcl[,3:5],1,which.max)
            v_gcl[,3:5]<-unmap(c_final2,group=1:L)
            
            
            if(s>1){
                change_measure[s] <- (mean(as.matrix(abs(u_gk[,2:4] - u_gk_old[,2:4]),na.rm=TRUE)))
                change_measure2[s] <- (mean(as.matrix(abs(v_gcl[,3:5] - v_gcl_old[,3:5]))))
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
        
        ### Getting the approximated estimates for u_gk and v_gcl
        
        combined_df <- do.call(rbind, u_gk_s)
        u_gk <- aggregate(. ~ Gene, data = combined_df, FUN = sum)
        u_gk<-u_gk[order(u_gk$Gene),]
        u_gk[,2:4]<-u_gk[,2:4]/length(u_gk_s)
        u_gk[,2:4]<-u_gk[,2:4]/rowSums(u_gk[,2:4])
        
        combined_df <- do.call(rbind, v_gcl_s)
        v_gcl <- aggregate(. ~ Gene+CpG, data = combined_df, FUN = sum)
        v_gcl<-v_gcl[order(v_gcl$Gene,v_gcl$CpG),]
        v_gcl[,3:5]<-v_gcl[,3:5]/length(v_gcl_s)
        v_gcl[,3:5]<-v_gcl[,3:5]/rowSums(v_gcl[,3:5])
        
        ## New dataframe with 
        u_gk_est<-u_gk
        v_gcl_est<-v_gcl
        
        ################## M-step #########################
        u_gk_est_old<-u_gk_est
        v_gcl_est_old<-v_gcl_est
        g_final2<-apply(u_gk[,2:4],1,which.max)
        u_gk_est[,2:4]<-unmap(g_final2,group=1:K)
        u_gk_est<-as.data.frame(u_gk_est)
        colnames(u_gk_est)<-c("Gene","1","2","3")
        cols_to_convert <- c("1", "2","3")
        u_gk_est[cols_to_convert] <- lapply(u_gk_est[cols_to_convert], as.numeric)
        c_final2<-apply(v_gcl[,3:5],1,which.max)
        v_gcl_est[,3:5]<-unmap(c_final2,group=1:L)
        v_gcl_est<-as.data.frame(v_gcl_est)
        colnames(v_gcl_est)<-c("Gene","CpG","1","2","3")
        cols_to_convert <- c("1", "2","3")
        v_gcl_est[cols_to_convert] <- lapply(v_gcl_est[cols_to_convert], as.numeric)
        
        
        tau_est<-colSums(u_gk_est[,2:(K+1)])/G
        
        pi_lk_est <- sapply(1:K, function(k) {
            sapply(1:L, function(l) {
                sum1 <- sum2 <- 0
                x=sapply(1:G, function(g) {
                    df_new <- v_gcl_est[v_gcl_est$Gene == u_gk_est$Gene[g], 3:(L+2)]
                    index <- which(v_gcl_est$Gene == u_gk_est$Gene[g])
                    Cg_num <- nrow(df_new)
                    if(Cg_num != 0)
                    {
                        sum1 <-  (u_gk_est[g, k + 1] * sum(df_new[, l]))
                        sum2 <-  (u_gk_est[g, k + 1] * Cg_num) 
                    }
                    c(sum1,sum2)
                })
                sum(x[1,])/sum(x[2,])
            })
        })
        
        
        mu_est<-vector()
        sigma_sq_est<-vector()
        for(k in 1:K){
            temp=0
            for(g in 1:G){
                temp_seq<-as.data.frame(seq_transformed_df[seq_transformed_df$Gene==u_gk_est$Gene[g], -1])
                for(n in 1:N){
                    temp<-temp+(u_gk_est[g, k + 1] * temp_seq[,n])
                }
            }
            mu_est[k]<-(temp)/(N*sum(u_gk_est[, k + 1]))
        }
        
        
        
        for(k in 1:K){
            temp=0
            for(g in 1:G){
                temp_seq<-as.data.frame(seq_transformed_df[seq_transformed_df$Gene==u_gk_est$Gene[g], -1])
                for(n in 1:N){
                    temp<-temp+(u_gk_est[g, k + 1] * ((temp_seq[,n] - mu_est[k])^2))
                }
            }
            sigma_sq_est[k]<-(temp)/(N*sum(u_gk_est[, k + 1]))
        }
        sigma_est<-sigma_sq_est^0.5
        
        sigma_temp<-sum(tau_est*sigma_est)
        sigma_est<-rep(sigma_temp,times=K)
        
        
        lambda_est <- sapply(1:L, function(l) {
            sum1 <- sum(sapply(1:C, function(c) {
                sum(v_gcl_est[c, l + 2] * meth_transformed_df[meth_transformed_df$CpG==v_gcl_est$CpG[c],3:(N+2)])
            }))
            sum1 / (N * sum(v_gcl_est[, l + 2]))
        })
        
        rho_sq_est <- numeric(L)
        
        for (l in 1:L) {
            Gc<-unique(v_gcl_start$Gene)
            sum2 <- numeric(length(Gc))
            sum4<-numeric(length(Gc))
            for (g in 1:length(Gc)) {
                cpg <- meth_transformed_df[meth_transformed_df$Gene== Gc[g], ]
                v_gcl_temp <- v_gcl_est[v_gcl_est$Gene == Gc[g], ]
                sum1 <- numeric(nrow(cpg))
                for (c in 1:nrow(cpg)) {
                    temp <- numeric(N)
                    for (n in 1:N) {
                        temp[n] <- v_gcl_temp[c, l + 2] * ((cpg[cpg$CpG==v_gcl_temp$CpG[c], n + 2] - lambda_est[l])^2)
                    }
                    sum1[c] <- sum(temp)
                }
                sum4[g] <- sum(v_gcl_temp[, l + 2])
                sum2[g] <- sum(sum1)
            }
            rho_sq_est[l] <- sum(sum2) / (N * sum(sum4))
        }
        
        rho_est<-rho_sq_est^0.5
        pi_temp<-rowSums(t(pi_lk_est)*(tau_est))
        rho_temp<-sum(pi_temp*rho_est)
        rho_est<-rep(rho_temp,times=L)
        
        v_gcl_est_old <- v_gcl_est_old[match(df_complete$CpG, v_gcl_est_old$CpG), ]
        class_gene=as.factor(apply(u_gk_est_old[,2:(K+1)], 1, which.max))
        class_cpg=as.factor(apply(v_gcl_est_old[,3:(L+2)], 1, which.max))
        
        
        ############ Convergence check #########################
        
        pi_diff<-mean(pi_lk_est-pi_lk)
        tau_diff<-mean(tau_est-tau)
        mu_diff<-mean(mu_est-mu)
        sigma_diff<-mean(sigma_est-sigma)
        lambda_diff<-mean(lambda_est-lambda)
        rho_diff<-mean(rho_est-rho)
        total_diff[i]<-(pi_diff+tau_diff+mu_diff+
                            sigma_diff+
                            lambda_diff +rho_diff)/6
        ########### Update ###########################
        pi_lk<-pi_lk_est
        tau<-tau_est
        mu<-mu_est
        sigma<-sigma_est
        lambda<-lambda_est
        rho<-rho_est
        
        if(abs(total_diff[i])<epsilon){
            break
        }
        
        
    }
    
    true_deg<-ifelse(df_gene$Group_Gene=="1",1,ifelse(df_gene$Group_Gene=="3",1,0))
    true_dmc<-ifelse(df_complete_raw$Group_CpG=="1",1,ifelse(df_complete_raw$Group_CpG=="3",1,0))
    ##########################Idiffomix##################################
    intmod_deg<-ifelse(class_gene=="1",1,ifelse(class_gene=="3",1,0))
    conf_matrix<-table(true_deg,intmod_deg)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_gene_intmod[seed,1]<-fdr
    perf_metrics_gene_intmod[seed,2]<-sensitivity
    perf_metrics_gene_intmod[seed,3]<-specificity
    perf_metrics_gene_intmod[seed,4]<-ari
    
    intmod_dmc<-ifelse(class_cpg=="1",1,ifelse(class_cpg=="3",1,0))
    conf_matrix<-table(true_dmc,intmod_dmc)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_cpg_intmod[seed,1]<-fdr
    perf_metrics_cpg_intmod[seed,2]<-sensitivity
    perf_metrics_cpg_intmod[seed,3]<-specificity
    perf_metrics_cpg_intmod[seed,4]<-ari

    #################### Mclust ##############################
    ## Mclust to get initial independent clustering
    gene_mc_1<-Mclust(seq_transformed_df[,2:5],K,modelNames = c("EVI","EII","EEI"))
    gene_cluster<-gene_mc_1$classification
    ls_index<-order(order(gene_mc_1$parameters$mean[1,]))
    gene_cluster_mclust<-ls_index[gene_cluster]
    cpg_mc_1<-Mclust(meth_transformed_df[,3:6],L,modelNames =c("EVI","EII","EEI"))
    cpg_cluster<-cpg_mc_1$classification
    ls_index<-order(order(cpg_mc_1$parameters$mean[1,]))
    cpg_cluster_mclust<-ls_index[cpg_cluster]
    mclust_deg<-ifelse(gene_cluster_mclust=="1",1,ifelse(gene_cluster_mclust=="3",1,0))
    
    conf_matrix<-table(true_deg,mclust_deg)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_gene_mclust[seed,1]<-fdr
    perf_metrics_gene_mclust[seed,2]<-sensitivity
    perf_metrics_gene_mclust[seed,3]<-specificity
    perf_metrics_gene_mclust[seed,4]<-ari

    mclust_dmc<-ifelse(cpg_cluster_mclust=="1",1,ifelse(cpg_cluster_mclust=="3",1,0))
    conf_matrix<-table(true_dmc,mclust_dmc)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_cpg_mclust[seed,1]<-fdr
    perf_metrics_cpg_mclust[seed,2]<-sensitivity
    perf_metrics_cpg_mclust[seed,3]<-specificity
    perf_metrics_cpg_mclust[seed,4]<-ari


  #####################Limma#############################
    limma_gene<-gene_cpm
    d0 <- DGEList(limma_gene)
    d0 <- calcNormFactors(d0)
    y <- voom(d0, design, plot = T)
    fit3 <- lmFit(y, design)
    fit4 <- eBayes(fit3)
    deff.gene <- topTable(fit4,adjust.method = "BH",
                          coef = 5,
                          number = Inf,sort.by = "none")
    d1_sorted <- deff.gene[match(rownames(seq_transformed_df), rownames(deff.gene)), ]
    d1_sorted <- d1_sorted[rownames(seq_transformed_df), ]
    
    gene_limma<-ifelse(d1_sorted$adj.P.Val<=0.05,1,0)
    
    fit <- lmFit(df_complete_raw[,c(3:10)], design)
    fit2 <- eBayes(fit)
    deff.meth = topTable(fit2, coef = 5,number = Inf, adjust.method = "BH",sort.by = "none")
    cpg_limma<-ifelse(deff.meth$adj.P.Val<=0.05,1,0)
    
    conf_matrix<-table(true_deg,gene_limma)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_gene_limma[seed,1]<-fdr
    perf_metrics_gene_limma[seed,2]<-sensitivity
    perf_metrics_gene_limma[seed,3]<-specificity
    perf_metrics_gene_limma[seed,4]<-ari

    conf_matrix<-table(true_dmc,cpg_limma)
    tn <- conf_matrix[1, 1]
    fp <- conf_matrix[1, 2]
    fn <- conf_matrix[2, 1]
    tp <- conf_matrix[2, 2]

    fdr <- fp / (fp + tp)
    specificity <-tn/(fp+tn)
    sensitivity <-tp/(tp+fn)
    ari<-unlist(classAgreement(conf_matrix)[4])

    perf_metrics_cpg_limma[seed,1]<-fdr
    perf_metrics_cpg_limma[seed,2]<-sensitivity
    perf_metrics_cpg_limma[seed,3]<-specificity
    perf_metrics_cpg_limma[seed,4]<-ari
    
    
    
    ####### ARI for case 3
    # conf_matrix<-table(intmod_deg,mclust_deg)
    # perf_metrics_gene_intmod[seed,1]<-unlist(classAgreement(conf_matrix)[4])
    # conf_matrix<-table(intmod_deg,gene_limma)
    # perf_metrics_gene_intmod[seed,2]<-unlist(classAgreement(conf_matrix)[4])
    # 
    # conf_matrix<-table(intmod_dmc,mclust_dmc)
    # perf_metrics_cpg_intmod[seed,1]<-unlist(classAgreement(conf_matrix)[4])
    # conf_matrix<-table(intmod_dmc,cpg_limma)
    # perf_metrics_cpg_intmod[seed,2]<-unlist(classAgreement(conf_matrix)[4])
    
}



