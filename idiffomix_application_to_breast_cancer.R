### install required packages
library(fourPNO)    
library(mclust)
library(readr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(limma)
library(missMethyl)
library(purrr)
library(e1071)
library(ggplot2)
library(edgeR)
library(MASS)
library(reshape2)
library(cowplot)
library(foreach)
library(doParallel)
library(circlize)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(grid)


###################### Data preprocessing ##########################
## Read data

load("gene_df.rda")    ## RNA-Seq data
load("meth_df.rda")    ## Methylation array data
load("legacy_df.rda")  ## Legacy data containing CHR,
                                                  ## Gene and MAPINFO

## Step - 1
## Select CpG sites having genes at the below positions along genome

relevant_positions <- c("TSS200", "5'UTR", "1stExon")
# Function to filter and select relevant genes
filter_genes <- function(gene_string, position_string) {
    genes <- unlist(strsplit(gene_string, ";"))
    positions <- unlist(strsplit(position_string, ";"))
    relevant_genes <- genes[positions %in% relevant_positions]
    if (length(relevant_genes) > 0) {
        return(relevant_genes[1])
    } else {
        return(NA)  # Return NA if no relevant genes are found
    }
}

# Apply the function and create a new dataframe with relevant genes

legacy_df$Relevant_genes <- mapply(filter_genes, legacy_df$UCSC_RefGene_Name, 
                                   legacy_df$UCSC_RefGene_Group)
legacy_df <- legacy_df[!is.na(legacy_df$Relevant_genes), ]



# Create the final dataframe

final_df <- legacy_df[, c("IlmnID","CHR","MAPINFO", "Relevant_genes")]
colnames(final_df)<-c("CpG","CHR","MAPINFO","Gene")
final_df=final_df[complete.cases(final_df$Gene),]
chr_gene<-as.vector(unique(final_df$Gene))



### Filter and transform Gene data

rownames(gene_df)<-gene_df$gene_id
gene_df_new<-gene_df[,-1]
targets <- data.frame(Patients=c(1,2,3,4,5,1,2,3,4,5),
                      Conditions=c("Benign","Benign",
                                   "Benign","Benign","Benign","Tumour","Tumour",
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
nonzero <- rowSums(dge$counts) > 5
dge %>% .[nonzero,]
dge %>% calcNormFactors
cpm_values <- cpm(dge)

# Transform filtered RNA-Seq count data

log_cpm_values <- log2(cpm_values + 1)  # Adding 1 to avoid log of zero
log_cpm_values<-as.data.frame(log_cpm_values)

## Step - 2

chr_cpg<-as.vector(unique(final_df$CpG))
meth_df_genome<-meth_df[meth_df$CpG %in% chr_cpg,]
meth_df_genome<-meth_df_genome[,c(1:11)]
meth_df_genome=merge(meth_df_genome,final_df,by.x = "CpG",by.y="CpG")
meth_df_genome<-meth_df_genome[,c(14,1,12,13,2:11)]
colnames(meth_df_genome)[2]<-"CpG"
meth_df_genome <- subset(meth_df_genome, CHR != "X" & CHR != "Y")

## Step - 3
## Filter CpG sites from Step - 2

log_cpm_values<-log_cpm_values[order(rownames(log_cpm_values)),]
meth_df_genome<-meth_df_genome[order(meth_df_genome$Gene,meth_df_genome$CHR,meth_df_genome$MAPINFO),]


## Calculating differential values for RNA-Seq data

seq_transformed_df<-log_cpm_values
seq_transformed_df<-as.data.frame(cbind(rownames(log_cpm_values),seq_transformed_df))
colnames(seq_transformed_df)[1]<-"Gene"
cols_to_convert <- c("Patient1_Cond1","Patient2_Cond1","Patient3_Cond1",
                      "Patient5_Cond1","Patient6_Cond1",
                     "Patient1_Cond2","Patient2_Cond2","Patient3_Cond2",
                     "Patient5_Cond2","Patient6_Cond2")
seq_transformed_df[cols_to_convert] <- lapply(seq_transformed_df[cols_to_convert], as.numeric)

seq_transformed_df[,2]<-(seq_transformed_df[,7]-seq_transformed_df[,2])
seq_transformed_df[,3]<-(seq_transformed_df[,8]-seq_transformed_df[,3])
seq_transformed_df[,4]<-(seq_transformed_df[,9]-seq_transformed_df[,4])
seq_transformed_df[,5]<-(seq_transformed_df[,10]-seq_transformed_df[,5])
seq_transformed_df[,6]<-(seq_transformed_df[,11]-seq_transformed_df[,6])
seq_transformed_df<-seq_transformed_df[,-c(7:11)]
seq_transformed_df<-as.data.frame(seq_transformed_df)
colnames(seq_transformed_df)<-c("Gene","Patient1","Patient2","Patient3","Patient5","Patient6")
cols_to_convert <- c("Patient1","Patient2","Patient3","Patient5","Patient6")
seq_transformed_df[cols_to_convert] <- lapply(seq_transformed_df[cols_to_convert], as.numeric)


## Calculating differential values for Methylation data

meth_df_genome<-meth_df_genome[meth_df_genome$Gene %in% seq_transformed_df$Gene,]
meth_transformed_df<-meth_df_genome
meth_transformed_df[,c(5:14)]<-log(meth_transformed_df[,5:14]/(1-meth_transformed_df[,5:14]))
log_meth<-meth_transformed_df
limma_meth<-meth_transformed_df
meth_transformed_df[,5]<-(meth_transformed_df[,10]-meth_transformed_df[,5])
meth_transformed_df[,6]<-(meth_transformed_df[,11]-meth_transformed_df[,6])
meth_transformed_df[,7]<-(meth_transformed_df[,12]-meth_transformed_df[,7])
meth_transformed_df[,8]<-(meth_transformed_df[,13]-meth_transformed_df[,8])
meth_transformed_df[,9]<-(meth_transformed_df[,14]-meth_transformed_df[,9])
meth_transformed_df<-meth_transformed_df[,-c(3,4,10:14)]

colnames(meth_transformed_df)<-c("Gene","CpG","Patient1","Patient2","Patient3","Patient5","Patient6")
cols_to_convert <- c("Patient1","Patient2","Patient3","Patient5","Patient6")
meth_transformed_df[cols_to_convert] <- lapply(meth_transformed_df[cols_to_convert], as.numeric)


load("gene_names.rda")
distinct_genes<-gene_names
distinct_genes<-distinct_genes[!(distinct_genes$CHR %in% c("X","Y")),]

################### Fitting joint mixture model with transformed data ##########################################

## set seed
set.seed(09082024)

## 
seq_transformed_df_genome<-seq_transformed_df
meth_transformed_df_genome<-meth_transformed_df

seq_transformed_df_genome<-seq_transformed_df_genome[seq_transformed_df_genome$Gene %in% distinct_genes$Gene,]
meth_transformed_df_genome<-meth_transformed_df_genome[meth_transformed_df_genome$Gene %in% distinct_genes$Gene,]
## Run the model for each chromosome individually and store results

chr_vec<-sort(as.numeric(unique(distinct_genes$CHR)))
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


################### Function initialization ####################
initialization<-function(seq_data,meth_data,K,L,N,probs=c(0.1,0.9)){
    gene_quantile<-quantile(as.vector(unlist(seq_data[,2:(N+1)])),probs = probs)
    gcluster_temp<-ifelse(rowMeans(seq_data[,2:(N+1)])<=gene_quantile[1],1,
                          ifelse(rowMeans(seq_data[,2:(N+1)])<gene_quantile[2],
                                 2,3))
    cpg_quantile<-quantile(as.vector(unlist(meth_data[,3:(N+2)])),probs = probs)
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
e_step<-function(seq_data,meth_data,u_gk,v_gcl,K,L,N,G,C,
                 tau,pi_lk,mu,sigma,lambda,rho,conv_th=0.001){
    
    # print(dim(pi_lk))
    ################ E-step ###################################
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
                        (dnorm(seq_data[g, n + 1], mu[k],
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
            temp3=apply(prob_u_gk[,2:4],1,max)
            prob_u_gk[,2:4]=exp(prob_u_gk[,2:4]-(temp3))
            prob_u_gk[,2:4]<-prob_u_gk[,2:4]/rowSums(prob_u_gk[,2:4])
            u_gk_s[[s]]<-prob_u_gk
            g_final2<-apply(prob_u_gk[,2:4],1,which.max)
            u_gk[,2:4]<-mclust::unmap(g_final2,group=1:K)
            
            
            # ####### v_gcl estimation ##############################
            prob_v_gcl<-matrix(0,nrow=C,ncol=L+2)
            for(c in 1:C){
                for(l in 1:L){
                    temp1 <- sum(sapply(1:N, function(n)
                        (dnorm(meth_data[c,n+2],lambda[l], 
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
            temp3=apply(prob_v_gcl[,3:5],1,max)
            prob_v_gcl[,3:5]=exp(prob_v_gcl[,3:5]-(temp3))
            prob_v_gcl[,3:5]<-prob_v_gcl[,3:5]/rowSums(prob_v_gcl[,3:5])
            v_gcl_s[[s]]<-prob_v_gcl
            c_final2<-apply(prob_v_gcl[,3:5],1,which.max)
            v_gcl[,3:5]<-mclust::unmap(c_final2,group=1:L)
            
            
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
        
        return(list(u_gk=u_gk,v_gcl=v_gcl))
}


######## Parallel code #############
# Register the number of cores you want to use
num_cores <- detectCores() - 1  # Adjust as needed
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


start_time=Sys.time()

# Parallelized foreach loop
results <- foreach(chr = chr_vec, .combine = 'list', .packages = c('mclust')) %dopar% {
    
    ## Select data for the corresponding CHR
    gene_chr <- as.vector(unlist(distinct_genes[distinct_genes$CHR == chr, "Gene"]))
    seq_transformed_df <- seq_transformed_df_genome[seq_transformed_df_genome$Gene %in% gene_chr, ]
    meth_transformed_df <- meth_transformed_df_genome[meth_transformed_df_genome$Gene %in% gene_chr, ]
    seq_transformed_df <- seq_transformed_df[order(seq_transformed_df$Gene), ]
    meth_transformed_df <- meth_transformed_df[order(meth_transformed_df$Gene, meth_transformed_df$CpG), ]
    
    ####################### Coordinate ascent method #########################
    ############## Initialization #######################
    K = 3
    L = 3
    N = 5
    G = nrow(seq_transformed_df)
    C = nrow(meth_transformed_df)
    iter = 20
    epsilon = 1e-5
    
    ## Quantile clustering
    initial_clustering <- initialization(seq_transformed_df, meth_transformed_df, K, L, N, probs = c(0.1, 0.9))
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
        g_final2 <- apply(u_gk[, 2:4], 1, which.max)
        u_gk[,2:4] <- as.data.frame(mclust::unmap(g_final2, group = 1:K))
        colnames(u_gk) <- c("Gene", "Cluster1", "Cluster2","Cluster3")
        
        c_final2 <- apply(v_gcl[, 3:5], 1, which.max)
        v_gcl[,3:5] <- as.data.frame(mclust::unmap(c_final2, group = 1:L))
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
    tau_final_par[[chr]] <- final_results[[i]]$tau_final
    mu_final_par[[chr]] <- final_results[[i]]$mu_final
    sigma_final_par[[chr]] <- final_results[[i]]$sigma_final
    pi_lk_final_par[[chr]] <- final_results[[i]]$pi_lk_final
    lambda_final_par[[chr]] <- final_results[[i]]$lambda_final
    rho_final_par[[chr]] <- final_results[[i]]$rho_final
    gene_cluster_final_par[[chr]] <- final_results[[i]]$gene_cluster_final
    cpg_cluster_final_par[[chr]] <- final_results[[i]]$cpg_cluster_final
    u_gk_final_par[[chr]] <- final_results[[i]]$u_gk_final
    v_gcl_final_par[[chr]] <- final_results[[i]]$v_gcl_final
}

end_time=Sys.time()

total_time=difftime(end_time,start_time,units = "mins")

final_gene_df_par<-as.data.frame(matrix(NA,nrow=0,ncol=ncol(seq_transformed_df)+1))
colnames(final_gene_df_par)<-c(colnames(seq_transformed_df),"Final_Cluster")
final_cpg_df_par<-as.data.frame(matrix(NA,nrow=0,ncol=ncol(meth_transformed_df)+1))
colnames(final_cpg_df_par)<-c(colnames(meth_transformed_df),"Final_Cluster")
for(chr in chr_vec){
    print(chr)
    # chr=22
    ## Select data for the corresponding CHR
    
    gene_chr<-as.vector(unlist(distinct_genes[distinct_genes$CHR==chr,"Gene"]))
    seq_transformed_df_2<-seq_transformed_df_genome[seq_transformed_df_genome$Gene %in% gene_chr,]
    meth_transformed_df2<-meth_transformed_df_genome[meth_transformed_df_genome$Gene %in% gene_chr,]
    seq_transformed_df_2<-seq_transformed_df_2[order(seq_transformed_df_2$Gene),]
    meth_transformed_df2<-meth_transformed_df2[order(meth_transformed_df2$Gene,
                                                     meth_transformed_df2$CpG),]
    seq_transformed_df_2$Final_Cluster<-gene_cluster_final_par[[chr]]
    meth_transformed_df2$Final_Cluster<-cpg_cluster_final_par[[chr]]
    final_gene_df_par<-rbind(final_gene_df_par,seq_transformed_df_2)
    final_cpg_df_par<-rbind(final_cpg_df_par,meth_transformed_df2)
    
}



########## Run mclust on whole genome
### Run mclust on whole together and do kegg and Go
seq_transformed_df_genome_2<-seq_transformed_df_genome
meth_transformed_df_genome_2<-meth_transformed_df_genome
seq_transformed_df_genome_2<-seq_transformed_df_genome_2[order(seq_transformed_df_genome_2$Gene),]
meth_transformed_df_genome_2<-meth_transformed_df_genome_2[order(meth_transformed_df_genome_2$Gene,meth_transformed_df_genome_2$CpG),]
mclust_gene_genome<-Mclust(seq_transformed_df_genome_2[,2:6],G=K,modelNames = c("EVI","EEI","EII"))
mclust_cpg_genome<-Mclust(meth_transformed_df_genome_2[,3:7],G=L,modelNames = c("EVI","EEI","EII"))

gene_cluster<-mclust_gene_genome$classification
ls_index<-order(order(mclust_gene_genome$parameters$mean[1,]))
mclust_gene_cluster<-ls_index[gene_cluster]
ls_z_index<-order(mclust_gene_genome$parameters$mean[1,])
z_gene_reordered <- mclust_gene_genome$z[, ls_z_index]

cpg_cluster<-mclust_cpg_genome$classification
ls_index<-order(order(mclust_cpg_genome$parameters$mean[1,]))
mclust_cpg_cluster<-ls_index[cpg_cluster]
ls_z_index<-order(mclust_cpg_genome$parameters$mean[1,])
z_cpg_reordered <- mclust_cpg_genome$z[, ls_z_index]


seq_transformed_df_genome_2<-as.data.frame(cbind(seq_transformed_df_genome_2,mclust_gene_cluster))
meth_transformed_df_genome_2<-as.data.frame(cbind(meth_transformed_df_genome_2,mclust_cpg_cluster))
colnames(seq_transformed_df_genome_2)[7]<-"Mclust_Cluster"
colnames(meth_transformed_df_genome_2)[8]<-"Mclust_Cluster"
final_gene_df_par<-merge(final_gene_df_par,seq_transformed_df_genome_2[,c(1,7)],by="Gene")
final_cpg_df_par<-merge(final_cpg_df_par,meth_transformed_df_genome_2[,c(1,2,8)],by=c("Gene","CpG"))

########## Run limma on whole genome

seq_limma_3<-gene_cpm[rownames(gene_cpm) %in% seq_transformed_df_genome$Gene,]
seq_limma_3<-seq_limma_3[order(rownames(seq_limma_3)),]
d0 <- DGEList(seq_limma_3)
d0 <- calcNormFactors(d0)
y <- voom(d0, design, plot = T)
fit3 <- lmFit(y, design)
fit4 <- eBayes(fit3)
deff.gene <- topTable(fit4,adjust.method = "BH",
                      coef = 6,
                      number = Inf,sort.by = "none")
d1_sorted <- deff.gene[match(final_gene_df_par$Gene, rownames(deff.gene)), ]
gene_limma<-ifelse(d1_sorted$adj.P.Val<=0.05,1,0)
gene_integromix<-ifelse(final_gene_df_par$Final_Cluster==1,1,ifelse(final_gene_df_par$Final_Cluster==3,1,0))
gene_mclust<-ifelse(final_gene_df_par$Mclust_Cluster==1,1,ifelse(final_gene_df_par$Mclust_Cluster==3,1,0))



meth_limma_3<-meth_df_genome[meth_df_genome$CpG %in% meth_transformed_df_genome$CpG,]
rownames(meth_limma_3)<-meth_limma_3$CpG
meth_limma_3<-meth_limma_3[,5:14]
meth_limma_3<-log(meth_limma_3/(1-meth_limma_3))
fit <- lmFit(meth_limma_3, design)
fit2 <- eBayes(fit)
deff.meth = topTable(fit2, coef = 6,number = Inf, adjust.method = "BH",sort.by = "none")
dcpg_sorted <- deff.meth[match(final_cpg_df2$CpG, rownames(deff.meth)), ]
cpg_limma<-ifelse(dcpg_sorted$adj.P.Val<=0.05,1,0)
cpg_integromix<-ifelse(final_cpg_df_par$Final_Cluster==1,1,ifelse(final_cpg_df_par$Final_Cluster==3,1,0))
cpg_mclust<-ifelse(final_cpg_df_par$Mclust_Cluster==1,1,ifelse(final_cpg_df_par$Mclust_Cluster==3,1,0))
final_gene_df_par$Limma_cluster<-gene_limma
final_cpg_df_par$Limma_cluster<-cpg_limma


##############################################################################
## Run Bayes' theorem fro graphs
tau= tau_final_par[[7]]
pi=pi_lk_final_par[[7]]
temp <- matrix(tau, 3, 3, byrow = TRUE) * pi
p_gene_given_cpg <- temp / rowSums(temp)
round(p_gene_given_cpg, 3)
P_gene_given_cpg <- p_gene_given_cpg


# Plotting the results
colnames(P_gene_given_cpg) <- c("E-", "E0", "E+")
rownames(P_gene_given_cpg) <- c("M-", "M0", "M+")
df_conditional <- melt(P_gene_given_cpg)
colnames(df_conditional) <- c("CpG_Cluster", "Gene_Cluster", "Probability")

# Plotting the results
df_prior <- data.frame(Cluster = as.factor( c("E-", "E0", "E+")), Probability = tau)
colnames(pi) <- c("E-", "E0", "E+")
rownames(pi) <- c("M-", "M0", "M+")
df_pi <- melt(pi)
colnames(df_pi) <- c("CpG_Cluster", "Gene_Cluster", "Probability")

df_prior$Cluster <- factor(df_prior$Cluster, levels = c("E-", "E0", "E+"))

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
final_plot <- grid.arrange(
    arrangeGrob(p1_labeled, p2_labeled, p3_labeled, nrow = 1),
    legend,
    ncol = 1,
    heights = c(10, 1)  # Adjust the height ratio if needed
)


######## Gene and meth graphs
#TNFRSF18 = 1, GPX7 = 1, RAD51 = 15
gene_name="RAD51"
chr=15
gene_reg_plot<-final_gene_df_par[final_gene_df_par$Gene==gene_name,]
cpg_reg_plot<-final_cpg_df_par[final_cpg_df_par$Gene==gene_name,]
cpg_order_plot<-meth_df_genome[meth_df_genome$Gene==gene_name,2:4]
merged_df <- merge(cpg_reg_plot, cpg_order_plot, by = "CpG")
ordered_df <- merged_df[order(merged_df$CHR, merged_df$MAPINFO), ]
ordered_df<-ordered_df[,c(2,1,3:10)]
cpg_reg_plot<-ordered_df

max_lim=max(max(gene_reg_plot[,2:7]),max(cpg_reg_plot[,3:7]))+0.2
min_lim=min(min(gene_reg_plot[,2:7]),min(cpg_reg_plot[,3:7]))-0.2

gene_reg_plot_long <- gene_reg_plot %>%
    pivot_longer(cols = starts_with("Patient"), names_to = "Patient", values_to = "Value") %>%
    mutate(Type = "Gene Expression")

# Prepare the CpG methylation data
cpg_reg_plot_long <- cpg_reg_plot %>%
    pivot_longer(cols = starts_with("Patient"), names_to = "Patient", values_to = "Value") %>%
    mutate(Type =  CpG)

# Combine the two datasets
combined_df <- bind_rows(gene_reg_plot_long, cpg_reg_plot_long)
combined_df$Final_Cluster<-ifelse(combined_df$Type=="Gene Expression","Gene",combined_df$Final_Cluster)
combined_df$Final_Cluster <- factor(combined_df$Final_Cluster, levels = c("Gene",1, 2,3), labels = c("Gene","M-", "M0","M+"))
gene_exp<-paste0(gene_name) 
combined_df$Type<-ifelse(combined_df$Type=="Gene Expression",
                         gene_exp,combined_df$Type)
combined_df$Type <- factor(combined_df$Type, levels = unique(combined_df$Type))
# Plotting
p1 <- ggplot(combined_df, aes(x = Type, y = Value, fill = Final_Cluster)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 1.5, color = "black") + # Add jittered points
    scale_fill_manual(values = c(
        "Gene"="gray",
        "M-" = "orange", "M0" = "#009E73", "M+" = "brown")) + # Set fill colors
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
p1

z_gene<-z_gene_reordered[rownames(z_gene_reordered)==gene_name,]


df_updated <- data.frame(Cluster = as.factor(c("E-","E0","E+")), Probability = z_gene)
df_updated$Cluster <- factor(df_updated$Cluster, levels = c("E-","E0","E+"))
p2 <- ggplot(df_updated, aes(x=Cluster, y=Probability,  fill = Probability)) +
    geom_bar(stat="identity", color="black", size=1.2)  +
    scale_fill_gradient(low="white", high="blue")+
    xlab("Gene cluster")+
    ylim(0, 1.1) +   theme_minimal()+theme(legend.position = "none",   # Ensure the plot background is clear
                                           axis.text.x = element_text(face = "bold", size = 16),  # Make x-axis ticks bold
                                           axis.text.y = element_text(face = "bold", size = 16), 
                                           axis.title.x = element_text(face = "bold", size = 16),  # Make x-axis label bold
                                           axis.title.y = element_text(face = "bold", size = 16)   )+
    geom_text(aes(label=sprintf("%.2f", Probability)), 
              vjust=-0.5, 
              size=6,              # Increase the text size
              fontface="bold")

u_gk_plot<-u_gk_final_par[[chr]]
p_gene_given_all<-unlist(u_gk_plot[u_gk_plot$Gene==gene_name,2:4])
df_posterior <- data.frame(Cluster = as.factor(c("E-","E0","E+")), Probability = p_gene_given_all)
df_posterior$Cluster <- factor(df_posterior$Cluster, levels = c("E-","E0","E+"))

p3 <- ggplot(df_posterior, aes(x=Cluster, y=Probability,  fill = Probability)) +
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
labelC <- textGrob("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))

# Arrange p1 and p2 with labels
p1_labeled <- arrangeGrob(p1, top = labelA)
p2_labeled <- arrangeGrob(p2, top = labelB)

# Combine p1 and p2 into p6
p3_labeled <- arrangeGrob(p3, top = labelC)
# p4_labeled <- arrangeGrob(p4, top = labelC)

# Arrange all plo#ts into a grid
final_plot <- grid.arrange(
    p1_labeled, 
    arrangeGrob(p2_labeled, p3_labeled, ncol = 2), 
    nrow = 2
)


##################### GO and KEGG analysis ##############
###### DMCs identified by idiffomix
#########################################
sigcpgs=final_cpg_df_par[final_cpg_df_par$Final_Cluster=="1"|
                          final_cpg_df_par$Final_Cluster=="3","CpG"]
allcpgs=final_cpg_df_par$CpG
gst_idiffomix <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO",
                         plot.bias = TRUE, prior.prob = TRUE)
table(gst_idiffomix$FDR<0.05)
topGSA(gst_idiffomix[gst_idiffomix$ONTOLOGY=="BP",],n=10)
gst_idiffomix_df<-gst_idiffomix[gst_idiffomix$FDR<0.05,]

kegg_idiffomix <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs,
                          collection = "KEGG",  sig.genes = TRUE)
topGSA(kegg_idiffomix, n=10)
table(kegg_idiffomix$FDR<0.05)
kegg_idiffomix_df=kegg_idiffomix[kegg_idiffomix$FDR<0.05,]



###### DMCs identified by mclust
##########################################
sigcpgs=final_cpg_df_par[final_cpg_df_par$Mclust_Cluster=="1"|
                          final_cpg_df_par$Mclust_Cluster=="3","CpG"]
allcpgs=final_cpg_df_par$CpG
gst_mclust <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO",
                     plot.bias = TRUE, prior.prob = TRUE)
table(gst_mclust$FDR<0.05)
topGSA(gst_mclust)
gst_mclust_df<-gst_mclust[gst_mclust$FDR<0.05,]

kegg_mclust <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs,
                      collection = "KEGG",  sig.genes = TRUE)
topGSA(kegg_mclust, n=5)
table(kegg_mclust$FDR<0.05)
kegg_mclust_df=kegg_mclust[kegg_mclust$FDR<0.05,]



###### DMCs identified by limma
##########################################
sigcpgs=final_cpg_df_par[cpg_limma==1,"CpG"]
allcpgs=final_cpg_df_par$CpG
gst_limma <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO",
                    plot.bias = TRUE, prior.prob = TRUE)
table(gst_limma$FDR<0.05)
topGSA(gst_limma)
gst_limma_df<-gst_limma[gst_limma$FDR<0.05,]

kegg_limma <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs,
                     collection = "KEGG",  sig.genes = TRUE)
topGSA(kegg_limma, n=5)
table(kegg_limma$FDR<0.05)
kegg_limma_df=kegg_limma[kegg_limma$FDR<0.05,]


########## DEGs by integromix #############
deg_genes <- final_gene_df_par[final_gene_df_par$Final_Cluster=="1" |
                                final_gene_df_par$Final_Cluster=="3","Gene"]

for(i in 1:length(unmapped_genes)){
    print(unmapped_genes[i])
    deg_genes[deg_genes==unmapped_genes[i]]=gene_replacement[gene_replacement$Old_genes==unmapped_genes[i],"New_genes"]
}
deg_entrez <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
unmapped_genes <- setdiff(deg_genes, deg_entrez$SYMBOL)
print(unmapped_genes)

# Perform GO enrichment analysis for degs by integrated model
go_enrichment <- enrichGO(gene          = deg_entrez$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

barplot(go_enrichment, showCategory = 10, title = "GO Enrichment Analysis")
head(go_enrichment)
go_enrichment_df <- as.data.frame(go_enrichment)


kegg_enrichment <- enrichKEGG(gene= deg_entrez$ENTREZID,
                              organism = 'hsa', # 'hsa' for human
                              pvalueCutoff = 0.05)

barplot(kegg_enrichment, showCategory = 10, title = "KEGG Enrichment Analysis")
head(kegg_enrichment)
kegg_enrichment_df <- as.data.frame(kegg_enrichment)

########## DEGs by mclust #############
deg_genes <- final_gene_df_par[final_gene_df_par$Mclust_Cluster=="1" |
                                final_gene_df_par$Mclust_Cluster=="3","Gene"]


for(i in 1:length(unmapped_genes)){
    print(unmapped_genes[i])
    deg_genes[deg_genes==unmapped_genes[i]]=gene_replacement[gene_replacement$Old_genes==unmapped_genes[i],"New_genes"]
}
deg_entrez <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
unmapped_genes <- setdiff(deg_genes, deg_entrez$SYMBOL)
print(unmapped_genes)


# Perform GO enrichment analysis for degs by integrated model
go_enrichment_mclust <- enrichGO(gene          = deg_entrez$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 keyType       = "ENTREZID",
                                 ont           = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05)

barplot(go_enrichment_mclust, showCategory = 10, title = "GO Enrichment Analysis")
head(go_enrichment_mclust)
go_enrichment_mclust_df <- as.data.frame(go_enrichment_mclust)


kegg_enrichment_mclust <- enrichKEGG(gene= deg_entrez$ENTREZID,
                                     organism = 'hsa', # 'hsa' for human
                                     pvalueCutoff = 0.05)

barplot(kegg_enrichment_mclust, showCategory = 10, title = "KEGG Enrichment Analysis")
head(kegg_enrichment_mclust)
kegg_enrichment_mclust_df <- as.data.frame(kegg_enrichment_mclust)

########## DEGs by limma #############
deg_genes <- final_gene_df_par[final_gene_df_par$Limma_cluster=="1","Gene"]

# gene_replacement[335,]<-c("TMEM49","VMP1")
for(i in 1:length(unmapped_genes)){
    print(unmapped_genes[i])
    deg_genes[deg_genes==unmapped_genes[i]]=gene_replacement[gene_replacement$Old_genes==unmapped_genes[i],"New_genes"]
}
deg_entrez <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
unmapped_genes <- setdiff(deg_genes, deg_entrez$SYMBOL)
print(unmapped_genes)


# Perform GO enrichment analysis for degs by integrated model
go_enrichment_limma <- enrichGO(gene          = deg_entrez$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "ENTREZID",
                                ont           = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05)

barplot(go_enrichment_limma, showCategory = 10, title = "GO Enrichment Analysis")
head(go_enrichment_limma)
go_enrichment_limma_df <- as.data.frame(go_enrichment_limma)

kegg_enrichment_limma <- enrichKEGG(gene= deg_entrez$ENTREZID,
                                    organism = 'hsa', # 'hsa' for human
                                    pvalueCutoff = 0.05)

barplot(kegg_enrichment_limma, showCategory = 10, title = "KEGG Enrichment Analysis")
head(kegg_enrichment_limma)
kegg_enrichment_limma_df <- as.data.frame(kegg_enrichment_limma)



