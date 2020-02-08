library("ggplot2")
library("devtools")
library("dplyr")
library("stringi")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# library("heatmap3")     #Load latest version of heatmap.3 function

sifrianm <- "D:\\HoxCancer\\Xena\\Samplist\\"
#sifrianm <- "C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\Xena\\Samplist\\"

files <- list.files(path=sifrianm, pattern="*.csv")
BRCA_filename <- "TCGA_GTEX_Breast_BRCA.csv"
print(files)

# out_file_type determines the type of the output file
#
out_file_type = 'pdf'
out_file_type = 'tiff'

linenum = 0
memuza_sum_mat = matrix(data=NA,ncol=5,nrow = 36)
colnames(memuza_sum_mat) = c("Cohort","CancerType","cophenetic mean","cophenetic sd","# of samples")


    #1  sifrianofile <- "D:\\HoxCancer\\Xena\\SamplistTest\\"
  sifria <- paste(sifrianm,BRCA_filename,sep="")
  
#  sifria <- paste("D:\\HoxCancer\\ZincFinger\\Samplist\\",filenm,sep="")
#  sifria <- paste("D:\\HoxCancer\\Ribosomal\\Samplist\\",filenm,sep="")
  cancertype_allcol <- read.table(sifria, header = TRUE, sep = ',',row.names=2)
  cancertype_allcol["HS"] <- stri_sub(cancertype_allcol$Samp.ID,-2,-1)
  TCGA_health <- cancertype_allcol[(cancertype_allcol$HS=="11" & cancertype_allcol$Cohort=="TCGA"),]
  TCGA_sick <- cancertype_allcol[(cancertype_allcol$HS=="01" & cancertype_allcol$Cohort=="TCGA"),]
  GTEX_health <- cancertype_allcol[(cancertype_allcol$Cohort=="GTEX"),]
  cancertype <- cancertype_allcol[,-1:-5]     # remove char columns
  cancertype_hox <- cancertype_allcol[,-1:-5]
  cancertype_hox <- cancertype[,-40:-63]

   cancertype <- as.matrix(cancertype)
  cancertype_hox <- as.matrix(cancertype_hox)

  cancertype_m<-apply(cancertype,c(1,2),as.numeric)
  cancertype_hox_m <- apply(cancertype_hox,c(1,2),as.numeric)

  #Set a working directory for output files
 
  #setwd("C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\R_cluster\\HoxCan")
  
  #Create vectors with column names and row names
  
  rown <- rownames(cancertype_m)             
  coln <- colnames(cancertype_m)
  
  ###Create color side bars
  cohort_v <- cancertype_allcol$Cohort
  rowClusters <- substr(cohort_v, 1, 1)
  rowCluster <- ifelse(rowClusters == 'G','lightgoldenrod1','green4')
  
  
  #Define custom dist and hclust functions for use with heatmaps
  mydist=function(c) {dist(c,method="euclidian")}
#  clust_method <- "median"
  clust_method <- "average"
  myclust=function(c) {hclust(c,method=clust_method)}

    #Create heatmap using custom heatmap.3 source code loaded above for all hoxes

 wdir <- paste(sifrianm,"BRCAHM\\",sep="")
 ##  setwd("D:\\HoxCancer\\TCGA\\HM")
 setwd(wdir)
##########################################################
  coln <- c(colnames(cancertype_hox_m))
  rownm <- c(rownames(cancertype_hox_m))
  
  gtex_rownmlist <-  sample(c(rownames(GTEX_health)),100,replace = FALSE)
  tcga_rownmlist <-   sample(c(rownames(TCGA_sick)),100,replace = FALSE)
  tcga_health_rownmlist <-   sample(c(rownames(TCGA_health)),100,replace = FALSE)
  
  tcga_hox_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% tcga_rownmlist, ] 
 
  tcga_hox_health_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% tcga_health_rownmlist, ]
  row.names(tcga_hox_health_m) <- paste("H",row.names(tcga_hox_health_m), sep = "")
  #   since there are less tcga samples need to sample same number from gtex
  gtex_hox_samp_m <- cancertype_hox_m[row.names(cancertype_hox_m) %in% gtex_rownmlist, ]
  gtex_sample_rownmlist <- row.names(gtex_hox_samp_m)
  # create a matrix that has same number of sick and same number of health samples
  TCGATCGA_hox_m <- rbind(tcga_hox_m,tcga_hox_health_m)
  TCGAGTEX_hox_m <- rbind(tcga_hox_m,gtex_hox_samp_m)
  
##########################################################  
  
  cancer_name <- substr(BRCA_filename,11,nchar(BRCA_filename)-4)
  
  if (out_file_type == 'pdf') {
    fn = paste("hm_TCGATCGA",cancer_name,".pdf",sep="")
    pdf(file=fn)
  } else 
    if (out_file_type == 'tiff') {
      fn = paste("hm_TCGATCGA",cancer_name,".tiff",sep="")
      tiff(filename = fn, height = 6.2, width = 6.2, units = 'in', res = 300)
    }
  
  #print(filenm)
#  second_title_line = paste(" HOX exp (method: ",clust_method,")")
#  second_title_line = paste(" RandRibos exp (method: ",clust_method,")")
  second_title_line <- paste("(TCGA:TCGA)" )
  underline_start_pos <- stri_locate_all(pattern = '_', cancer_name, fixed = TRUE)
  last_underline_pos <- underline_start_pos[[1]][length(underline_start_pos[[1]])]
  cancer_name_acronym <- substr(cancer_name,last_underline_pos+1,stri_length(cancer_name))
  first_title_line= paste("HOX exp", cancer_name_acronym, "- Healthy vs. Tumor")
  first_title_line= paste("HOX exp", cancer_name_acronym, "-")
  second_title_line <- paste("Healthy(TCGA) vs. Tumor(TCGA)" )
  
  main_title= c(first_title_line,' ',second_title_line)
  par(cex.main=1)

  colnSnH <- c(colnames(TCGATCGA_hox_m))
  
  # for (o_rn in tcga_health_rownmlist)
  # {
  #   n_rn <- paste('H',o_rn,sep="")
  #   names(TCGATCGA_hox_m)[which(names(TCGATCGA_hox_m)==o_rn)]= n_rn
  # }
  rownmSnH <- c(rownames(TCGATCGA_hox_m))
#  row.names(x) <- value
  
  rowClusters <- substr(rownmSnH, 1, 1)
  rowCluster <- ifelse(rowClusters == 'H','lightgoldenrod1','green4')
  legendTitle<- ifelse(rowClusters == 'H','Healthy (GTEX)','Tumor (TCGA)') # to be used for the legend 
  
  row_annotation <- as.matrix(t(rowCluster))
  rownames(row_annotation) <- c("H/T")
  can_category <- as.factor

  heatmap.3(TCGATCGA_hox_m, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
            Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, labRow=TRUE, cexCol = 0.5,cexRow=0.05,
            density.info="none", trace="none", main=main_title, labCol=colnames(TCGATCGA_hox_m), col=rev(heat.colors(75)),
            ColSideColorsSize=7,RowSideColors=row_annotation, RowSideColorsSize=1, 
            KeyValueName="exp lvl")
 # legend("topright", legend = c("Tumor (TCGA)","Healthy (TCGA)"), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  legend("topright", legend = unique(legendTitle), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  dev.off()
#**************************************************************************************************
  if (out_file_type == 'pdf') {
    fn = paste("hm_TCGAGTEX",cancer_name,".pdf",sep="")
    pdf(file=fn)
  } else 
    if (out_file_type == 'tiff') {
      fn = paste("hm_TCGAGTEX",cancer_name,".tiff",sep="")
      tiff(filename = fn, height = 6.2, width = 6.2, units = 'in', res = 300)
    }
  #print(filenm)
#  second_title_line = paste(" HOX exp (method: ",clust_method,")")
  #  second_title_line = paste(" RandRibos exp (method: ",clust_method,")")
  second_title_line <- paste("(GTEX:TCGA)" )
  first_title_line= paste("HOX exp", cancer_name_acronym, "- Healthy vs. Tumor")
  
  first_title_line= paste("HOX exp", cancer_name_acronym, "-")
  second_title_line <- paste("Healthy(GTEX) vs. Tumor(TCGA)" )
  
  main_title= c(first_title_line,' ',second_title_line)
  par(cex.main=1)
  
  colnSnH <- c(colnames(TCGAGTEX_hox_m))
  rownmSnH <- c(rownames(TCGAGTEX_hox_m))
  
  
  rowClusters <- substr(rownmSnH, 1, 1)
  rowCluster <- ifelse(rowClusters == 'G','lightgoldenrod1','green4')
  legendTitle<- ifelse(rowClusters == 'G','Healthy (GTEX)','Tumor (TCGA)') # to be used for the legend
  
  row_annotation <- as.matrix(t(rowCluster))
  rownames(row_annotation) <- c("H/T")
  can_category <- as.factor
  
  heatmap.3(TCGAGTEX_hox_m, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
            Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, labRow=TRUE, cexCol = 0.5,cexRow=0.05,
            density.info="none", trace="none", main=main_title, labCol=colnames(TCGAGTEX_hox_m), col=rev(heat.colors(75)),
            ColSideColorsSize=7,RowSideColors=row_annotation, RowSideColorsSize=1, 
            KeyValueName="exp lvl")
 # legend("topright", legend = c("Tumor (TCGA)","Healthy (GTEX)"), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  legend("topright", legend = unique(legendTitle), col = unique(rowCluster),  lty= 1, lwd = 5, cex=.5)
  
  dev.off()
##}
