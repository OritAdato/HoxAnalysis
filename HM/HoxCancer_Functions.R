## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   
## Fuction that counts and returns the number of color changes in the heatmap side color bar
##
number_of_ColChanges<- function(hc,tcgaonly){
  num_of_colChange =0
  first_sampInOrder = hc[["order"]] [1]
  prev_samp_name = hc[["labels"]] [first_sampInOrder]
  if (tcgaonly==TRUE){
    prev_samp_first = hc[["labels"]] [first_sampInOrder]
    prev_samp_first_let = substr(hc[["labels"]] [first_sampInOrder],nchar(prev_samp_first)-1,
                                 nchar(prev_samp_first)-1)
  }  else{
    prev_samp_first_let = substr(hc[["labels"]] [first_sampInOrder],1,1)
  }
  
  for (sample_num in  hc[["order"]]) 
  {
   
    curr_samp_name = hc[["labels"]] [sample_num]
    if (tcgaonly==TRUE){
      curr_samp_first_let <- substr(curr_samp_name,nchar(curr_samp_name)-1,
                                   nchar(curr_samp_name)-1)
    }  else{
      curr_samp_first_let <- substr(hc[["labels"]] [sample_num],1,1)
    }
    
    stam1 <- c('prev:', prev_samp_name, 'curr:',curr_samp_name)
#    print (stam1)
    if (curr_samp_first_let!= prev_samp_first_let)
    {
      num_of_colChange = num_of_colChange + 1
      stam <- c(curr_samp_name,prev_samp_first_let,curr_samp_first_let,sample_num,num_of_colChange)
      prev_samp_first_let = curr_samp_first_let
      prev_samp_name = curr_samp_name
#      print(stam)
    }
  }
  print(num_of_colChange)
  return(num_of_colChange)
}

## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  
calc_avg_dist <- function(cancertype_hox_m_SnH,submat_rn){
  
  healthy_p <- cancertype_hox_m_SnH[submat_rn,]
  #m_testmatdist <- as.matrix(testmatdist)
  healthy_pdist <- mydist(healthy_p)
  m_healthy_pdist <- as.matrix(healthy_pdist)
  #colMeans(m_healthy_pdist)
  avg_dist <- mean(m_healthy_pdist)

  return(avg_dist)
}

## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
calc_avg_dist_ofHealthFromTumor <- function(cancertype_hox_m_SnH,GTEX_rn,TCGA_rn){
  
  vec_avg <- c()

  for (GTEX_samp in GTEX_rn) {
    for (TCGA_samp in TCGA_rn) {
      pair_samp <- rbind(cancertype_hox_m_SnH[GTEX_samp,],cancertype_hox_m_SnH[TCGA_samp,])
      pair_dist <- mydist(pair_samp)
      pair_dist_m <- as.matrix(pair_dist)
      vec_avg <- c(vec_avg, pair_dist_m[2,1])   # since it is distance between 2 vectors only the matrix is 2x2
    }
    
  }
  avg_dist_HealthfromTumor <- mean(vec_avg)
  return(avg_dist_HealthfromTumor)
}

######################################################################################################################################
# The function below compares the euclidan distance between healthy smaples to healthy samples and compares it with the avg eauclidean
# distance between healthy smaples to tumor samples 
######################################################################################################################################
ttest_avg_dist_ofHealthFromTumor <- function(cancertype_hox_m_SnH,GTEX_rn,TCGA_rn){
  
  vec_avgDist_HfromH <- c()   # to calulate average distand between healthy samples and themselves
  
  for (GTEX_samp in GTEX_rn) {
    for (GTEX_samp2 in GTEX_rn) {
      if (GTEX_samp == GTEX_samp2) {next()}    # not comparing sample with itself
      pair_samp <- rbind(cancertype_hox_m_SnH[GTEX_samp,],cancertype_hox_m_SnH[GTEX_samp2,])
      pair_dist <- mydist(pair_samp)
      pair_dist_m <- as.matrix(pair_dist)
      vec_avgDist_HfromH <- c(vec_avgDist_HfromH, pair_dist_m[2,1])   # since it is distance between 2 vectors only the matrix is 2x2
    }
  }
  vec_avgDist_HfromT <- c()  # to calulate average distand between healthy samples and tumor samples
  
  for (GTEX_samp in GTEX_rn) {
    for (GTEX_TCGA in TCGA_rn) {
      #   if (GTEX_samp == GTEX_samp2) {next()}
      pair_samp <- rbind(cancertype_hox_m_SnH[GTEX_samp,],cancertype_hox_m_SnH[GTEX_TCGA,])
      pair_dist <- mydist(pair_samp)
      pair_dist_m <- as.matrix(pair_dist)
      vec_avgDist_HfromT <- c(vec_avgDist_HfromT, pair_dist_m[2,1])   # since it is distance between 2 vectors only the matrix is 2x2
    }
  }
  
  tt_out <- t.test(vec_avgDist_HfromH,vec_avgDist_HfromT)
  
  avg_dist_HealthfromTumor <- mean(vec_avgDist_HfromT)
  avg_dist_HealthfromHealth <- mean(vec_avgDist_HfromH)
  ci_h <- tt_out$conf.int[1]
  ci_ht <-tt_out$conf.int[2]
  outp_l <- list(tt_out$p.value,tt_out$estimate,ci_h,ci_ht)
  
  return(outp_l)
}
