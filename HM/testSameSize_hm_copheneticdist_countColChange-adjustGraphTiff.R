#library("gplots")
library("ggplot2")
library("devtools")
library("dplyr")
library("stringi")
#library("heatmap3")
#source("D:\\HoxCancer\\R_cluster\\HoxCancer_Functions.R")
source("C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\R_cluster\\HoxCancer_Functions.R")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#
# change: 19/12/2019
# create HM based on csv files and count number of side color changes the HM created in tiff format to change
# resulotion - need to change the 'res' parameter in the tiff command
#
#direprefix <- "D:\\"
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
#setwd("D:\\HoxCancer\\R_cluster\\CanHM\\Median")
#setwd("D:\\HoxCancer\\R_cluster\\CanHM")
linenum = 0
ttest_linenum = 0

memuza_sum_mat = matrix(data=NA,ncol=6,nrow = 60)
colnames(memuza_sum_mat) = c("Cohort","CancerType","cophenetic mean","cophenetic sd","# of samples",
                             "# of color changes")

# The matrix below will contain the result of the ttest performed for every ttest that is run
# for every gene comparing expression of tumor and health samples
#
ttest_results_mat = (matrix(data=NA,ncol = 6,nrow = 1000))
colnames(ttest_results_mat) = c("Cancer_type","Gene_name","ttest-p.val","mean_of_GTEX","mean_of_TCGA",
                                "Is_syg")

for(filenm in files)
 {
##  filenm <- 'TCGA_GTEX_Colon.csv'
  sifrianofile <- sifrianm
#4  sifrianofile <- "D:\\HoxCancer\\lsRibosomal\\Samplist\\"
#3  sifrianofile <- "D:\\HoxCancer\\Xena\\SamplistFPKM\\"
#2  sifrianofile <- "D:\\HoxCancer\\Xena\\Samplist\\"
#1  sifrianofile <- "D:\\HoxCancer\\Xena\\SamplistTest\\"

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
  #cancertype_mt <- t(cancertype_m)
  
  #Load latest version of heatmap.3 function
   
  #Set a working directory for output files
 
  #setwd("C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\R_cluster\\HoxCan")
  
  #Create vectors with column names and row names
  
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
  
 wdir <- paste(sifrianm,"HM19122019\\",sep="")
 setwd(wdir)
##########################################################
  coln <- c(colnames(cancertype_hox_m))
  rownm <- c(rownames(cancertype_hox_m))
  
  gtex_rownmlist <-  rownm[substr(rownm, 1, 4)=="GTEX"]
  tcga_rownmlist <-   rownm[substr(rownm, 1, 4)=="TCGA"]
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
  fn = paste("hm_",cancer_name,".tiff",sep="")
  #pdf(file=fn)
  tiff(filename = fn, height = 6.2, width = 6.2, units = 'in', res = 300)
  print(filenm)
  second_title_line = paste(gene_group,"exp (hclust method:",clust_method,")")
#  second_title_line = paste(" RandRibos exp (method: ",clust_method,")")
  main_title= c(cancer_name,second_title_line)
  par(cex.main=1)

  colnSnH <- c(colnames(cancertype_hox_m_SnH))
  rownmSnH <- c(rownames(cancertype_hox_m_SnH))
  
  ###Create color side bars  
  rowClusters <- substr(rownmSnH, 1, 1)
  rowCluster <- ifelse(rowClusters == 'G','lightgoldenrod1','green4')
  legendTitle<- ifelse(rowClusters == 'G','Healthy','Tumor') # to be used for the legend 
  
  row_annotation <- as.matrix(t(rowCluster))
  rownames(row_annotation) <- c("H/T")
  can_category <- as.factor
 # lmat = rbind(c(0,4,5),c(3,1,2))
#  lwid = c(1.5,4,1)
#  lhei = c(1.5,4)
  heatmap.3(cancertype_hox_m_SnH, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
            Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, labRow=TRUE, cexCol = 0.5,cexRow=0.05,
            density.info="none", trace="none", main=main_title, labCol=colnames(cancertype_m), col=rev(heat.colors(75)),
            ColSideColorsSize=3,RowSideColors=row_annotation, RowSideColorsSize=1, #lmat=lmat, lwid = lwid, lhei=lhei, 
            KeyValueName="Normalized exp lvl")
##  legend("topright", legend = unique(rowClusters), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  legend("topright", legend = unique(legendTitle), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  #d hm3_p = heatmap.3(cancertype_hox_m_SnH, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
  #d                   Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, labRow=TRUE, cexCol = 0.5,cexRow=0.05,
  #d                   density.info="none", trace="none", main=main_title, labCol=colnames(cancertype_m), col=rev(heat.colors(75)),
  #d                 ColSideColorsSize=3,RowSideColors=row_annotation, RowSideColorsSize=1, 
  #d                   KeyValueName="Normalized exp lvl")
  #d print(c(rownmSnH[rev(hm3_p$rowInd)]))   # print clustered rowid order
  #d print (c(colnSnH[hm3_p$colInd]))        # print clustered colid order
  dev.off()
  
  ######################################################################
  # For every gebe check if  the mean expression level of health samples is same as the mean
  # expression level of tumor samples
  
  GTEX_rn <- rownmSnH[substr(rownmSnH,1,4) == "GTEX"] 
  TCGA_rn <- rownmSnH[substr(rownmSnH,1,4) == "TCGA"]
  for (hox_gnm in colnSnH)
  {
    ttest_result<- t.test(cancertype_hox_m_SnH[GTEX_rn,hox_gnm],cancertype_hox_m_SnH[TCGA_rn,hox_gnm],paired=FALSE) 
    mean_of_x <- ttest_result$estimate[[1]]   # [[]] to removed the name from named value
    mean_of_y <- ttest_result$estimate[[2]]
    if (is.nan(ttest_result$p.value) == TRUE) 
    {
      is_syg <- "NaN" 
    } else if (ttest_result$p.value < 0.05) {
      is_syg <- "LT0.05" 
    } else {
      is_syg <- "NO"}
    
    ttest_result_line <- c(cancer_name, hox_gnm, ttest_result$p.value,mean_of_x[1],mean_of_y[1],is_syg)
    ttest_linenum <- ttest_linenum + 1
    ttest_results_mat[ttest_linenum,1:6] <- ttest_result_line[1:6]
  }
  #  tr<- t.test(cancertype_hox_m_SnH[GTEX_rn,"HOXD3"],cancertype_hox_m_SnH[TCGA_rn,"HOXD3"],paired=FALSE)
  #
  #############################################################################
  
#
# create the cophenetic matrix for cancer expression matrix that includes same number of 
# sick and health samples
#  
  hc <- hclust(mydist(cancertype_hox_m_SnH),method = clust_method)
##  print(hc$labels[c(hc$order)]) print order of labels
  
  linenum <- linenum + 1
  cop_hc <- as.matrix(cophenetic(hc))
  tcgaonly = FALSE
  n_of_ColChange <- number_of_ColChanges(hc,tcgaonly)
  tot_line <- c("GTEXnTCGA",cancer_name,mean(cop_hc[lower.tri(cop_hc)]),sd(cop_hc[lower.tri(cop_hc)]),dim(cancertype_hox_m_SnH)[1],
                n_of_ColChange)
  memuza_sum_mat[linenum,1:6] <- tot_line[1:6]
  
  linenum <- linenum + 1
  cop_hc_h <- cop_hc[substr(row.names(cop_hc),1, 4)=="GTEX",substr(colnames(cop_hc),1, 4)=="GTEX"]
  tot_line <- c("GTEX",cancer_name,mean(cop_hc_h[lower.tri(cop_hc_h)]),sd(cop_hc_h[lower.tri(cop_hc_h)]),dim(cancertype_hox_m_SnH)[1]/2)
  memuza_sum_mat[linenum,1:5] <- tot_line[1:5]
  
  linenum <- linenum + 1
  cop_hc_s <- cop_hc[substr(row.names(cop_hc),1, 4)=="TCGA",substr(colnames(cop_hc),1, 4)=="TCGA"]
  tot_line <- c("TCGA ",cancer_name,mean(cop_hc_s[lower.tri(cop_hc_s)]),sd(cop_hc_s[lower.tri(cop_hc_s)]),dim(cancertype_hox_m_SnH)[1]/2)
  memuza_sum_mat[linenum,1:5] <- tot_line[1:5]

 
}  
  dir_path_parts <- strsplit(shem_dir,"\\\\")
  ggroup <- dir_path_parts[[1]][length(dir_path_parts[[1]])-1]
  memuza_fn <- paste("memuza_sum_mat_",ggroup,".csv",sep = "")
  write.csv(memuza_sum_mat, file = memuza_fn, row.names = FALSE)
  ttest_fn <- paste("ttest_results_mat_",ggroup,".csv",sep = "")
  write.csv(ttest_results_mat, file = ttest_fn, row.names = FALSE)
  
#  write.csv(memuza_sum_mat, file = "memuza_sum_mat.csv", row.names = FALSE)
#  write.csv(ttest_results_mat, file = "ttest_results_mat.csv", row.names = FALSE)
  
  