#
# based on the program testSameSize_hm_copheneticdist_countColChange-adjustGraph.R
#
library("ggplot2")
library("devtools")
library("dplyr")
library("stringi")

#source("D:\\HoxCancer\\R_cluster\\HoxCancer_Functions.R")
source("C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\R_cluster\\HoxCancer_Functions.R")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#
# change: 02/10
# create HM based on csv files and count number of side color changes
#
direprefix <- "D:\\"
direprefix <- "C:\\Users\\orita\\Documents\\BIU\\"
##0411 shem_dir <- "HoxCancer\\ZincFinger\\Samplist\\"
####C:\Users\orita\Documents\BIU\HoxCancer\lsRibosomal\Samplist
shem_dir <- "HoxCancer\\Xena\\RANDgenes2\\Samplist\\"
shem_dir <- "HoxCancer\\Xena\\Samplist\\"     # for hox genes
#shem_dir <- "HoxCancer\\Xena\\SamplistTest\\"     # for hox genes
#shem_dir <- "HoxCancer\\lsRibosomal\\Samplist\\"     # for lsribosomal genes

gene_group = "HOX"
#gene_group = "Ribosomal"

sifrianm <- paste(direprefix,shem_dir,sep="")
#09/10 sifrianm <- paste(direprefix,"HoxCancer\\Ribosomal\\Samplist\\",sep="")
#09/10 sifrianm <- paste(direprefix,"HoxCancer\\Xena\\Samplist\\",sep="")

#4 files <- list.files(path=sifria, pattern="*.csv")
files <- list.files(path=sifrianm, pattern="*.csv")
print(files)

linenum = 0
ttest_linenum = 0
 
avg_euclidean_dist_mat = matrix(data=NA,ncol=403,nrow=14)
numOfRuns<- 100
###colnames(avg_euclidean_dist_mat) = c("Cohort","avg_GTEXdist","avg_TCGAdist","avg_dist_H_from_T")
coltitles <- c("Cohort")

# The matrix below will contain the result of the ttest performed for every ttest that is run
# for every gene comparing expression of tumor and health samples
#

filenum <- 0
for(filenm in files)
  {
  # filenm <- 'TCGA_GTEX_Colon_COAD.csv'
  # filenm <- 'TCGA_GTEX_brain_LGG.csv'
  filenum <- filenum + 1
  sifrianofile <- sifrianm

  sifria <- paste(sifrianofile,filenm,sep="")

#  sifria <- paste("D:\\HoxCancer\\ZincFinger\\Samplist\\",filenm,sep="")
#  sifria <- paste("D:\\HoxCancer\\Ribosomal\\Samplist\\",filenm,sep="")
  cancertype_allcol <- read.table(sifria, header = TRUE, sep = ',',row.names=2)
  cancertype_allcol <- cancertype_allcol[cancertype_allcol$Category != "No Category",]  # to remove TCGA health samples
  cancertype <- cancertype_allcol[,-1:-5]     # remove char columns
  cancertype_hox <- cancertype_allcol[,-1:-5]
  cancertype_hox <- cancertype[,-40:-62]
 

   cancertype <- as.matrix(cancertype)
  cancertype_hox <- as.matrix(cancertype_hox)

  cancertype_m<-apply(cancertype,c(1,2),as.numeric)
  cancertype_hox_m <- apply(cancertype_hox,c(1,2),as.numeric)

  rown <- rownames(cancertype_m) 
  coln <- colnames(cancertype_m)
  
  ###Create color side bars
 # cohort_v <- cancertype_allcol$Cohort
 # rowClusters <- substr(cohort_v, 1, 1)
 # rowCluster <- ifelse(rowClusters == 'G','lightgoldenrod1','green4')
  
  
  #Define custom dist and hclust functions for use with heatmaps
  mydist=function(c) {dist(c,method="euclidian")}
#  clust_method <- "median"
  clust_method <- "average"
  myclust=function(c) {hclust(c,method=clust_method)}

  # Create heatmap using custom heatmap.3 source code loaded above for all hoxes
  
 wdir <- paste(sifrianm,"HMeDist\\",sep="")
 setwd(wdir)
##########################################################
  coln <- c(colnames(cancertype_hox_m))
  rownm <- c(rownames(cancertype_hox_m))
  
  gtex_rownmlist <-  rownm[substr(rownm, 1, 4)=="GTEX"]
  tcga_rownmlist <-   rownm[substr(rownm, 1, 4)=="TCGA"]
  
  cancer_name <- substr(filenm,11,nchar(filenm)-4)
  avg_euclidean_dist_mat[filenum,1] <- c(cancer_name)

  GtexSetnum <- 2
  NextGtexSetnum <- 5
####-------> update number of runs  
  for (runn in 1:numOfRuns) {
# 100 samples from every type of cancer
  gtex_rownmlist <-  sample(rownm[substr(rownm, 1, 4)=="GTEX"],100,replace = FALSE)
  tcga_rownmlist <-   sample(rownm[substr(rownm, 1, 4)=="TCGA"],100,replace = FALSE)
  
 
  # cancertype_hox_m_sick <- cancertype_hox_m[row.names(cancertype_hox_m) %in% gtex_rownmlist,]
  # Below line builds a subset of distance matrix that includes only health lines

  if (length(gtex_rownmlist) < length(tcga_rownmlist)) {
    # number of samples in matrix is determined by the number of gtex sample which is the lower
    gtex_hox_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% gtex_rownmlist, ] 
    #   since there are less gtex samples need to sample same number from tcga 
    tcga_hox_m_samplnm <- sample(tcga_rownmlist,length(gtex_rownmlist), replace=FALSE)
    tcga_hox_samp_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% tcga_hox_m_samplnm, ]
    tcga_sample_rownmlist <- row.names(tcga_hox_samp_m)
    # create a matrix that has same number of sick and same number of health samples
    cancertype_hox_m_SnH <- rbind(tcga_hox_samp_m,gtex_hox_m)
    
  } else {
    # number of samples in cancer matrix is determined by the number of tcga samples which is the lower
    tcga_hox_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% tcga_rownmlist, ] 
    #   since there are less tcga samples need to sample same number from gtex 
    gtex_hox_m_samplnm <- sample(gtex_rownmlist,length(tcga_rownmlist), replace=FALSE)
    gtex_hox_samp_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% gtex_hox_m_samplnm, ]
    gtex_sample_rownmlist <- row.names(gtex_hox_samp_m)
    # create a matrix that has same number of sick and same number of health samples
    cancertype_hox_m_SnH <- rbind(gtex_hox_samp_m,tcga_hox_m)
    
  }

##########################################################  
  
  cancer_name <- substr(filenm,11,nchar(filenm)-4)
   ###print(filenm)

  colnSnH <- c(colnames(cancertype_hox_m_SnH))
  rownmSnH <- c(rownames(cancertype_hox_m_SnH))

  
  ######################################################################
  # For every gebe check if  the mean expression level of health samples is same as the mean
  # expression level of tumor samples
  
  GTEX_rn <- rownmSnH[substr(rownmSnH,1,4) == "GTEX"] 
  TCGA_rn <- rownmSnH[substr(rownmSnH,1,4) == "TCGA"]
 
# check average euclidean distance between tumor samples and 
# check average euclidean distance between healthy smaples
   
  avg_GTEX  <- calc_avg_dist(cancertype_hox_m_SnH,GTEX_rn)
  avg_TCGA  <- calc_avg_dist(cancertype_hox_m_SnH,TCGA_rn)
#  check average euclidean distance of every healthy samples from all other tumor samples
  
 ### avg_dist_ofHealthFromTumor <- calc_avg_dist_ofHealthFromTumor(cancertype_hox_m_SnH,GTEX_rn,TCGA_rn)
  comp_evg_euclideanDist_HfromT <- ttest_avg_dist_ofHealthFromTumor(cancertype_hox_m_SnH,GTEX_rn,TCGA_rn)
  tt_pval <- unname(comp_evg_euclideanDist_HfromT[[1]])
  avg_dist_HfromH <- unname(comp_evg_euclideanDist_HfromT[[2]][1])
  avg_dist_HfromT <- unname(comp_evg_euclideanDist_HfromT[[2]][2])
  ci_distHfromH <- unname(comp_evg_euclideanDist_HfromT[[3]])
  ci_distHfromT <- unname(comp_evg_euclideanDist_HfromT[[4]])
  #avg_euclidean_dist_mat[filenum,GtexSetnum:NextGtexSetnum] <- c(avg_GTEX,avg_dist_HfromH,avg_TCGA,avg_dist_HfromT,tt_pval,ci_distHfromH,ci_distHfromT)
  
  avg_euclidean_dist_mat[filenum,GtexSetnum:NextGtexSetnum] <- c(avg_TCGA,avg_GTEX,avg_dist_HfromT,tt_pval)
  
  coltitles <- append(coltitles,paste('R', runn, '_avg_TCGA',sep=""))
  coltitles <- append(coltitles,paste('R', runn, '_avg_GTEX',sep=""))
  #coltitles <- append(coltitles,paste('R', runn, '_avg_GTEX-F',sep=""))

  coltitles <- append(coltitles,paste('R', runn, '_avg_DistHfT',sep=""))
  coltitles <- append(coltitles,paste('R', runn, '_tt_pval',sep=""))
  #coltitles <- append(coltitles,paste('R', runn, '_ci_distHfromH',sep=""))
  #coltitles <- append(coltitles,paste('R', runn, '_ci_distHfromT',sep=""))
  
  #GtexSetnum <- GtexSetnum + 3
  #NextGtexSetnum <- NextGtexSetnum + 3
  
  GtexSetnum <- GtexSetnum + 4
  NextGtexSetnum <- NextGtexSetnum + 4
  
  ##print (c("avgGTEX",avg_GTEX,"avgTCGA",avg_TCGA,"avgHfromT",avg_dist_ofHealthFromTumor))
 
   #############################################################################
  
  }
  
}  

avg_euclidean_dist_df <- as.data.frame(avg_euclidean_dist_mat)
colnames(avg_euclidean_dist_df)[1:GtexSetnum-1] <- c(coltitles[1:GtexSetnum-1])

vec_col <- c()
for (rn in row.names(avg_euclidean_dist_df)){
  dfm <- avg_euclidean_dist_df %>% select(ends_with("_avg_GTEX"))
  dfm2 <- data.frame(sapply(dfm, function(x) as.numeric(as.character(x))))
  vec_col <- append(vec_col,mean(as.numeric(c(dfm2[rn,]))))
}
avg_euclidean_dist_df[[GtexSetnum]] <- vec_col
colnames(avg_euclidean_dist_df)[GtexSetnum] <- 'Mean_allruns_avg_GTEX'

vec_col <- c()
for (rn in row.names(avg_euclidean_dist_df)){
  dfm <- avg_euclidean_dist_df %>% select(ends_with("avg_DistHfT"))
  dfm3 <- data.frame(sapply(dfm, function(x) as.numeric(as.character(x))))
  vec_col <- append(vec_col,mean(as.numeric(c(dfm3[rn,]))))
}
avg_euclidean_dist_df[[GtexSetnum+1]] <- vec_col
colnames(avg_euclidean_dist_df)[GtexSetnum+1] <- 'Mean_allruns_avg_DistHfT'

  avg_eDist_fn <- paste("avg_Euclidean_distRuns",".csv",sep = "")
  write.csv(avg_euclidean_dist_df, file = avg_eDist_fn, row.names = FALSE)


  # dfm <- avg_euclidean_dist_df %>% select(ends_with("_avg_GTEX-F"))
  # dfm2 <- data.frame(sapply(dfm, function(x) as.numeric(as.character(x))))
  # vec <- as.numeric(c(dfm2[1,]))
  # mean(vec)
  # 