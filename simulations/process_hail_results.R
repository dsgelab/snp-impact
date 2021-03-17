library(data.table)

list_final <- NULL
RES_DESC <- NULL
for (i in 0:71)
{
  system(paste0("gsutil cp gs://mattia/output_sakari_sim/rsid_AF_beta_gold_pheno", i, "_* /Users/andreaganna/Downloads/"))
  system(paste0("gsutil cp gs://mattia/output_sakari_sim/rsid_beta_se_gwas_pheno", i, "* /Users/andreaganna/Downloads/"))
  
  file_o <- list.files(path="/Users/andreaganna/Downloads/","rsid_AF_beta_gold_pheno*")
  
  system(paste0("mv /Users/andreaganna/Downloads/", file_o, " /Users/andreaganna/Downloads/", gsub(".bgz",".gz",file_o)))
  system(paste0("mv /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,".tsv.bgz", " /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,".tsv.gz"))
  system(paste0("mv /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,"_logistic.tsv.bgz", " /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,"_logistic.tsv.gz"))
  
  file_o <- list.files(path="/Users/andreaganna/Downloads/","rsid_AF_beta_gold_pheno*")
  orig <- fread(paste0("/Users/andreaganna/Downloads/", file_o))
  new <- fread(paste0("/Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,".tsv.gz"))
  new_logistic <- fread(paste0("/Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,"_logistic.tsv.gz"))
  
  h2 <- as.numeric(gsub("h2","",grep("h2",strsplit(file_o,"_")[[1]], value=TRUE)))
  pi <- as.numeric(gsub("pi","",grep("pi",strsplit(file_o,"_")[[1]], value=TRUE)))
  k <- as.numeric(gsub("k|\\.tsv.gz","",grep("k",strsplit(file_o,"_")[[1]], value=TRUE)))
  
  df <- data.frame(rsid=orig$rsid,beta_gold=orig$beta,beta_gwas_linear=new$beta,se_gwas_linear=new$standard_error,pvalue_gwas_linear=new$p_value,beta_gwas_logistic=new_logistic$beta,se_gwas_logistic=new_logistic$standard_error,pvalue_gwas_logistic=new_logistic$p_value)
  
  list_final[[i+1]] <- df
  RES_DESC <- rbind(RES_DESC,c(i+1,h2,pi,k))
  system(paste0("rm /Users/andreaganna/Downloads/", file_o))
  system(paste0("rm /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,".tsv.gz"))
  system(paste0("rm /Users/andreaganna/Downloads/rsid_beta_se_gwas_pheno",i,"_logistic.tsv.gz"))
}


  beta_gold = do.call(cbind,sapply(list_final, "[", "beta_gold" ))
  pval_linear = do.call(cbind,sapply(list_final, "[", "pvalue_gwas_linear" ))
  pval_logistic = do.call(cbind,sapply(list_final, "[", "pvalue_gwas_logistic" ))
  
  # Select variants trurly associated with between 1 and 6 traits
  count_true_beta = which(rowSums(beta_gold!=0)<6 & rowSums(beta_gold!=0)>0)
  
  # Find which P-value is GW-significant
  count_GW_significant <- which(rowSums(pval_logistic < 0.00000005)>0)
  
  # Variants that are GW-significant in at least on trait and that the beta is trurly associated between 1 and 6 traits
  final_variants <- intersect(count_true_beta,count_GW_significant)
  
  # Select only the subset of variants
  beta_gold_subset <- sapply(list_final, function(x){x[final_variants,"beta_gold"]})
  pval_logistic_subset <- sapply(list_final, function(x){x[final_variants,"pvalue_gwas_logistic"]})
  
  ## Test if the overall the variants true betas have low p-values. That is the case
  for (i in 1:nrow(ff))
  {
    beta_gold_temp <- beta_gold_subset[i,]
    pval_logistic_temp <- pval_logistic_subset[i,]
    print(pval_logistic_temp[which(beta_gold_temp!=0)])
  }

  
