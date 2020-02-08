# package to write to xls: openxlsx
library("devtools")
library("survival")
library("survminer")
library("gtools")
library("xlsx")
library("stringi")

#
# Prepare list of gene pair of which their pair KM was more significant then the KM for every one of the genes separately
# The program runs for ++ ONE ++ cancaer type - file name is updated in the variable filenm
# The ouput xls sheet includes expression values of differnt medians (therefore program is called extended)
# To create the KM plots need to uncommet the section at the bottome below "Create KM plots" comment
#
# survival tutorial from: http://www.sthda.com/english/wiki/survival-analysis-basics
direprefix <- "D:\\HoxCancer\\Xena\\SamplistWsurvival\\"
#direprefix <- "C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\Xena\\SamplistWsurvival\\"
file_type = 'tiff'
file_type = 'pdf'
#inputdir <- "C:\\Users\\orita\\Documents\\BIU\\HoxCancer\\Xena\\SamplistWsurvival\\TestOut\\"
inputdir <- paste(direprefix,"TestOut\\",sep="")
###outdir <- paste(direprefix,"OutPairInOneGraph\\",sep="")
if (file_type == 'pdf'){
  outdir <- paste(direprefix,"PairsTrioGraph-PDFbonf\\",sep="")
} else {
  outdir <- paste(direprefix,"PairsTrioGraph-TiffBonf\\",sep="")
}
#outpltdir <- paste(outdir,"Pancreas-tiff\\",sep="")
#allSygKMpairsName_fn <- "testallsygKMAllcan1308AllinOne.xlsx"
allSygKMpairsName_fn <- "testallsygKMAllcan1308AllinOneWilcox.xlsx"     #12/01/2020


filenm <- 'TCGA_GTEX_Liver_LIHC_wSurvival.csv'
filenm <-'TCGA_GTEX_brain_GBM__wSurvival.csv'
filenm <-'TCGA_GTEX_brain_LGG__wSurvival.csv'
#filenm <- 'TCGA_GTEX_Leukemia_LAML_wSurvival.csv'
#filenm <- 'TCGA_GTEX_Stomach_STAD_wSurvival.csv'
#filenm <- 'TCGA_GTEX_Colon_COAD_wSurvival.csv'
#filenm <- 'TCGA_GTEX_Liver_LIHC_wSurvival.csv'
#filenm <- 'TCGA_GTEX_Lung_LUSC_wSurvival.csv'
#filenm <- 'TCGA_GTEX_Lung_LUAD_wSurvival.csv'
filenm <- 'TCGA_GTEX_Pancreas_PAAD_wSurvival.csv'

files <- list.files(path=direprefix, pattern="*.csv")
sifrianf_nm <- paste(direprefix,filenm,sep="")
cancertype_surv_df <- read.table(sifrianf_nm, header = TRUE, sep = ',',row.names=1)
colnm <- colnames(cancertype_surv_df)
colnm_HOX <- colnm[1:39]

setwd(inputdir)
KMpairsList_df<- read.xlsx(file=allSygKMpairsName_fn,sheetName = 'allSygKMallcan')
#setwd(outdir)
# 
# G1 = "HOXA1"
# G2 = "HOXA13"

df_allsyg_allCan <- data.frame(matrix(ncol = 45, nrow = 0))
df_allsygColTitle <- c("canType","G1", "G2", "Pair", "pvalNotTextBoth", "pvalBoth",
                       "hG1_maxBoth","hG1_minBoth","hG1_medB","tG1_maxBoth","tG1_minBoth","tG1_medB",
                       "hG2_maxBoth","hG2_minBoth","hG2_medB","tG2_maxBoth","tG2_minBoth","tG2_medB",
                       "numOf_onlyG1_HTMed","numOf_pnlyG1_LTMed","numof_BothG1G2_HTMed","numof_BothG1G2_LTMed",
                       "G1onlyKM_pval","G1only_pvalNum","Max_G1Only","Min_G1Only","Median_G1Only",
                       "G2only_pvalNum","G2onlyKM_pval","Max_G2Only","Min_G2Only","Median_G2Only",
                       "hG1p_max","hG1p_min","hG1p_med","tG1p_max","tG1p_min","tG1p_med",
                       "hG2p_max","hG2p_min","hG2p_med","tG2p_max","tG2p_min","tG2p_med","DoesBothChange")


colnames(df_allsyg_allCan) <- df_allsygColTitle

df_allsyg <- data.frame(matrix(ncol = 45, nrow = 0))
colnames(df_allsyg) <- df_allsygColTitle


p_HOX<- combinations(n=39,r=2,v=colnm_HOX,repeats.allowed=F)
p_HOX_rn = nrow(p_HOX)
bonferoni_corect_pval = 0.05/(39*38/2)

get_orderedDF_of_oneGene <- function(G1,cancertype_surv_df,len_bothHTM,len_bothLTM)
{
  ordered_by_G1 <- cancertype_surv_df[order(cancertype_surv_df[[G1]],decreasing = TRUE),]
  ordered_by_G1$gExp <- paste(substr(G1,4,nchar(G1)),"_expHTM",sep="")
  ordered_by_G1[(rownames(tail(ordered_by_G1,n=len_bothLTM))),"gExp"] <- paste(substr(G1,4,nchar(G1)),"_expLTM",sep="")
  G1_df <- rbind(head(ordered_by_G1,n=len_bothHTM), tail(ordered_by_G1,n=len_bothLTM))
  return(G1_df)
}
plot_graf <- function(fit_result,gtitle,legend_HTM,legend_LTM, pval2display){
  
  if (grepl("LGG" , gtitle)) {
    if ( stri_length(legend_HTM) > 20){
      legend_leftc <- 0.38
    } else {
      legend_leftc <- 0.24
    }
    legend_rightc <- 0.15
  } else {
    legend_leftc <- 0.66
    legend_rightc <- 0.97
    
  }
  
  p11 <-ggsurvplot(fit_result,
                    pval = pval2display, pval.size = 6, pval.coord = c(0, 0.03), conf.int = TRUE,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "strata", # Change line type by groups
                    title = gtitle, font.title = c(16, "bold", "darkblue"),
                    #        surv.median.line = "hv", # Specify median survival
                  #  legend = c(0.66, 0.97),  # determines the place of legend on the graph
                    legend = c(legend_leftc, legend_rightc),  # determines the place of legend on the graph
                   
                    legend.title = "geneExp",
                    legend.labs = c(legend_HTM,legend_LTM),
                    ggtheme = theme_classic(), # Change ggplot2 theme
                    break.time.by = 250, # # break X axis in time intervals by 250.
                    xlim = c(0,1780),    # # present narrower X axis
                    palette = c("red","navyblue"))
  p1 <- ggpar(p11, 
        font.title = c(16, "bold.italic", "black"),
        font.subtitle = c(14,  "blue"),
        font.caption = c(14,  "blue"),
        font.x = c(14,  "black"),
        font.y = c(14,  "black"),
        font.legend = 12,
        font.tickslab = 14
  )
  return(p1)
}


for (filenm in files){
sifrianf_nm <- paste(direprefix,filenm,sep="")
cancertype_surv_df <- read.table(sifrianf_nm, header = TRUE, sep = ',',row.names=1)

#cantype = substr(filenm,11,nchar(filenm)-14)
cantypeLongName = substr(filenm,11,nchar(filenm)-14)
if (grepl("brain" , filenm)){
  cantype = substr(filenm,11,nchar(filenm)-18)
} else {
  cantype = substr(filenm,11,nchar(filenm)-19)
}
underline_start_pos <- stri_locate_all(pattern = '_', cantypeLongName, fixed = TRUE)
last_underline_pos <- underline_start_pos[[1]][length(underline_start_pos[[1]])]
cancer_name_acronym <- substr(cantypeLongName,last_underline_pos+1,stri_length(cantypeLongName))

# 
# if (cantype == "brain" | cantype =="Breast" | cantype =="Colon" | cantype == "Esophagus") {
#   next()
# }


outpltdir <- paste(outdir,cantypeLongName,"\\",sep="")


for (j in 1:p_HOX_rn){
  G1 = p_HOX[j,1]
  G2 = p_HOX[j,2]
  

  a <- paste("--------->",G1,G2,cancer_name_acronym)
  print("*************************")
  print(a)
  print("*************************")
  
  both_HTM = paste(substr(G1,4,nchar(G1)),substr(G2,4,nchar(G2)),"_HTMed",sep="")
  both_LTM = paste(substr(G1,4,nchar(G1)),substr(G2,4,nchar(G2)),"_LTMed",sep="")
  
  cancertype_surv_df$bgExp <- paste(substr(G1,4,nchar(G1)),"_HTMed",sep="")
  cancertype_surv_df[cancertype_surv_df[[G1]] <= median(cancertype_surv_df[[G1]]),"bgExp"] <- 
      paste(substr(G1,4,nchar(G1)),"_LTMed", sep="")         # Now high median of G1 indicated as G1_HTMed and
                                                             # Low median of G1 indicated as G1_LTMed
 
  cancertype_surv_df$bgExp[(cancertype_surv_df[[G2]] <= median(cancertype_surv_df[[G2]]) & 
                              cancertype_surv_df[[G1]] <= median(cancertype_surv_df[[G1]]))] <- both_LTM

  cancertype_surv_df$bgExp[(cancertype_surv_df[[G2]] > median(cancertype_surv_df[[G2]]) & 
                                cancertype_surv_df[[G1]] > median(cancertype_surv_df[[G1]]))] <- both_HTM 
  
  len_bothHTM = length(cancertype_surv_df$bgExp[cancertype_surv_df$bgExp == both_HTM])
  len_bothLTM = length(cancertype_surv_df$bgExp[cancertype_surv_df$bgExp == both_LTM])
  len_G1HTM = length(which(cancertype_surv_df[["bgExp"]]== (paste(substr(G1,4,nchar(G1)),"_HTMed",sep=""))))
  len_G1LTM = length(which(cancertype_surv_df[["bgExp"]]== (paste(substr(G1,4,nchar(G1)),"_LTMed",sep=""))))
                           
                           
  if (len_bothLTM > 10 & len_bothHTM > 10) { 
    # print("*************************")
    # a <- paste("--------->",G1,G2)
    # print(a)
    both_LnH_df = subset(cancertype_surv_df,bgExp == both_LTM | bgExp == both_HTM,select = HOXA1:bgExp)
 
    fit <- survfit(Surv(OS.time, OS) ~ bgExp, data = cancertype_surv_df)
    fitboth <- survfit(Surv(OS.time, OS) ~ bgExp, data = both_LnH_df)
    # print(fit)
    # summary(fit)
    sv_fit <- surv_pvalue(fit)
    sv_fitboth <- surv_pvalue(fitboth)
    if (sv_fitboth["pval"] < bonferoni_corect_pval) {
      
      df1 <- get_orderedDF_of_oneGene(G1, cancertype_surv_df,len_bothHTM,len_bothLTM)
      df2 <- get_orderedDF_of_oneGene(G2, cancertype_surv_df,len_bothHTM,len_bothLTM)
      
      f_df1 <- survfit(Surv(OS.time, OS) ~ gExp, data = df1)
      f_df2 <- survfit(Surv(OS.time, OS) ~ gExp, data = df2)
      
      rishon <- surv_pvalue(f_df1)
      sheni <- surv_pvalue(f_df2)
      # 
      # rishon
      # sheni
      # sv_fitboth["pval"]
      
      if (sv_fitboth["pval"]< rishon["pval"] & sv_fitboth["pval"] < sheni["pval"]) { 
        print("both pval lower yes")
        print(rishon)
        print(sheni)
        print(sv_fitboth)
        a <- paste("--------->",G1,G2)
        print("*************************")
        print(a)
        print("$$$$   fit all    $$$$$$")
        print(fit[["strata"]])
        print (sv_fit)
        print("$$$$$$$$  fit both_HL  $$$$$$")
        print(fitboth[["strata"]])
        print(sv_fitboth)
        
        hdf1 <- head(df1,len_bothHTM) 
        hdf2 <- head(df2,len_bothHTM)
        tdf1 <- tail(df1,len_bothLTM)
        tdf2 <- tail(df2,len_bothLTM)
        
        max_hg1_only <- max(hdf1[[G1]])
        min_hg1_only <- min(hdf1[[G1]])
        max_tg1_only <- max(tdf1[[G1]])
        min_tg1_only <- min(tdf1[[G1]])
        
        max_hg2_only <- max(hdf2[[G2]])  
        min_hg2_only <- min(hdf2[[G2]])
        max_tg2_only <- max(tdf2[[G2]])
        min_tg2_only <- min(tdf2[[G2]])

        max_both_hg1 <- max(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_HTM])  #
        min_both_hg1 <- min(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_HTM])  #
        max_both_tg1 <- max(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_LTM])  #
        min_both_tg1 <- min(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_LTM])  #
        
        max_both_hg2 <- max(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_HTM])  #
        min_both_hg2 <- min(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_HTM])  #
        max_both_tg2 <- max(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_LTM])  #
        min_both_tg2 <- min(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_LTM])  #
        
        
        median_hg1_only <- median(hdf1[[G1]])
        median_tg1_only <- median(tdf1[[G1]])
        median_hg2_only <- median(hdf2[[G2]])
        median_tg2_only <- median(tdf2[[G2]])
        
 
        median_both_hg1 <- median(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_HTM]) #
        median_both_hg2 <- median(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_HTM]) #
        median_both_tg1 <- median(cancertype_surv_df[[G1]][cancertype_surv_df$bgExp == both_LTM]) #
        median_both_tg2 <- median(cancertype_surv_df[[G2]][cancertype_surv_df$bgExp == both_LTM]) #
        
         
        median_all_g1 <- median(cancertype_surv_df[[G1]])  #
        median_all_g2 <- median(cancertype_surv_df[[G2]])  #
        max_all_g1 <- max(cancertype_surv_df[[G1]])        #
        max_all_g2 <- max(cancertype_surv_df[[G2]])        #
        min_all_g1 <- min(cancertype_surv_df[[G1]])        #
        min_all_g2 <- min(cancertype_surv_df[[G2]])        #
 
  
    # building the xls sheet that include all details about the data     

        
        if ((grepl("brain" , filenm)) | (grepl("Lung" , filenm)) ){
          can_pair_name <- paste(cantypeLongName,G1,G2)
        } else {
          can_pair_name <- paste(cantype,G1,G2)
        }     
        doesboth_pairmemChange <- KMpairsList_df$PairChange[KMpairsList_df$CancNpair== can_pair_name]
        
        df_syg <- data.frame(
          canType = cancer_name_acronym,
          G1 = G1,
          G2 = G2,
          pair = paste(G1, G2),
          pvalNotTextBoth = unname(sv_fitboth["pval"]),
          pvalBoth = unname(sv_fitboth["pval.txt"]),
          hG1_maxBoth <- max_both_hg1,
          hG1_minBoth <- min_both_hg1,
          hG1_medB    <- median_both_hg1,
          tG1_maxBoth <- max_both_tg1,
          tG1_minBoth <- min_both_tg1,
          tG1_medB    <- median_both_tg1,
          
          hG2_maxBoth <- max_both_hg2,
          hG2_minBoth <- min_both_hg2,
          hG2_medB    <- median_both_hg2,
          tG2_maxBoth <- max_both_tg2,
          tG2_minBoth <- min_both_tg2,
          tG2_medB    <- median_both_tg2,
          

          numOf_onlyG1_HTMed= len_G1HTM,
          numOf_pnlyG1_LTMed= len_G1LTM,
          numof_BothG1G2_HTMed= len_bothHTM,
          numof_BothG1G2_LTMed= len_bothLTM,
          G1only_pval = unname(rishon["pval.txt"]),
          G1only_pvalNum = unname(rishon["pval"]),
          max_g1 <- max_all_g1,
          min_g1 <- min_all_g1,
          median_g1 <- median_all_g1,
  

          G2only_pvalNum = unname(sheni["pval"]),
          G2only_pval = unname(sheni["pval.txt"]),
          max_g2 <- max_all_g2,
          min_g2 <- min_all_g2,
          median_g2 <- median_all_g2,
          # hdf1 includes the head of the (high median of G1 ordered by expression level).Number of selected samples is the same as
          # the number of samples that appear in high median of both G1 and G2
          # similarly tdf1 includes the bottom median
          # and hdf2 includes same as above for G2
          #
          hG1p_max <- max_hg1_only,     
          hG1p_min <- min_hg1_only,
          hG1p_med <- median_hg1_only,
          tG1p_max <- max_tg1_only,
          tG1p_min <- min_tg1_only,
          tG1p_med <- median_tg1_only,
          
          hG2p_max <- max_hg2_only,     
          hG2p_min <- min_hg2_only,
          hG2p_med <- median_hg2_only,
          tG2p_max <- max_tg2_only,
          tG2p_min <- min_tg2_only,
          tG2p_med <- median_tg2_only,
          DoesBothChange <- doesboth_pairmemChange
          
        )
###=======================================================================================
        #   Create KM plots trios of pair and the separate genes
  
#     #   if (KMpairsList_df$DoesPairChangeTogether[KMpairsList_df$CancNpair== can_pair_name] == can_pair_name){
        if (KMpairsList_df$PairChange[KMpairsList_df$CancNpair== can_pair_name] == 'BothChange'){
 
          pval2display <- paste("pval=",toString(format.pval(sv_fitboth["pval"],digits = 3)),sep="")
          g_title <- paste(cancer_name_acronym, '-', G1,'&', G2)
          if (file_type == 'tiff'){
            fn_2saveplots <- paste(G1,'_',G2,'_',cancer_name_acronym,'.tiff',sep = "")     #for tiff
          } else {
            fn_2saveplots <- paste(G1,'_',G2,'_',cancer_name_acronym,'.pdf',sep = "")     #for pdf
          }
         #
          legend_HTM = paste(G1,'&',G2,' > medianExp',sep = "")
          legend_LTM = paste(G1,'&',G2,' < medianExp',sep = "")
          plt_both <- plot_graf(fitboth,g_title,legend_HTM,legend_LTM,pval2display)
          
          pval2display <- paste("pval=", toString(format.pval(rishon["pval"],digits = 3)),sep="")
          g_title <- paste(cancer_name_acronym, '-', G1)
          legend_HTM = paste(G1,' > medianExp',sep = "")
          legend_LTM = paste(G1,' < medianExp',sep = "")
          plt_G1 <- plot_graf(f_df1,g_title,legend_HTM,legend_LTM,pval2display)
          
          pval2display <-  paste("pval=",toString(format.pval(sheni["pval"],digits = 3)),sep="")
          g_title <- paste(cancer_name_acronym, '-', G2)
          legend_HTM = paste(G2,' > mediamExp',sep = "")
          legend_LTM = paste(G2,' < mediamExp',sep = "")
          plt_G2 <- plot_graf(f_df2,g_title,legend_HTM,legend_LTM,pval2display)
          plt_list <- list(plt_both, plt_G1, plt_G2)
         
        #   #arrange_ggsurvplots(plt_list, print = TRUE, ncol=3, risk.table.height = 0.15)
        #   # Arrange and save into pdf file
          setwd(outpltdir)
          res_plotlist <- arrange_ggsurvplots(plt_list, print = FALSE, ncol=3, risk.table.height = 0.15)
          if (file_type == 'tiff'){
            tiff(filename = fn_2saveplots, width =55, height = 30, units = "cm", res = 300)   #for tiff
            print(res_plotlist)                                                               #for tiff
 #         ggsave(fn_2saveplots, res_plotlist,width =55, height = 30, units = "cm")         #for pdf
            dev.off() #for tiff 
          } else {
            ggsave(fn_2saveplots, res_plotlist,width =55, height = 30, units = "cm")
          }
      }
###=======================================================================================
        df_allsyg <- rbind(df_allsyg,df_syg)
      }
    }
    
    }
}
}
 
 setwd(outdir)
 colnames(df_allsyg) <- df_allsygColTitle   # otherwhise titles include ....
# 
# sname <- paste("allSyg13_",cantype,sep="")
write.xlsx(df_allsyg, file = "testallsygKMpairs14012020.xlsx", 
            sheetName="allsygcan", append=TRUE)
# fit <- survfit(Surv(OS.time, OS) ~ bgExp, data = cancertype_surv_df)
# fitboth <- survfit(Surv(OS.time, OS) ~ bgExp, data = both_LnH_df)
# # print(fit)
# summary(fit)
# surv_pvalue(fit)
# surv_pvalue(fitboth)
##summary(fit)$table
# The function survfit() returns a list of variables, including the following components:
# The components can be accessed as follow:
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# d <- data.frame(time = fit$time,
#                 n.risk = fit$n.risk,
#                 n.event = fit$n.event,
#                 n.censor = fit$n.censor,
#                 surv = fit$surv,
#                 upper = fit$upper,
#                 lower = fit$lower
# )
# # head(d)
# # Change color, linetype by strata, risk.table color by strata
# ggsurvplot(fit,
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            break.time.by = 250, # # break X axis in time intervals by 250.
#            xlim = c(0,1780),    # # present narrower X axis
#            palette = c("#c00000", "#000080", "#006400", "#FF4500"))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# check values : 
# length(cancertype_surv_df$bgExp[cancertype_surv_df$bgExp == "OneHigherThenMedian"])
# length(cancertype_surv_df$bgExp[cancertype_surv_df$bgExp == "LowerThenMedian"])
# length(cancertype_surv_df$bgExp[cancertype_surv_df$bgExp == "bothlowerThenMedian"])
# check what values included in column: 
# levels(factor(cancertype_surv_df$bgExp))