from hoxcancer_functions import *
import gzip
import os.path
import random as rd
import pandas as pd
import json
from scipy import stats as sc
import numpy as np
from datetime import datetime
import pprint
#
# This program compares the expression level of random sets of 39 genes of Breast tissue -  TCGA:TCGA vs TCGA:GTEX
# Number of sets is controled by the variable: num_of_sets_in_oneRead
# To run this program for HOX set only need to uncomment : dict_of_genel_4run =...
#
def load_json2list(gene_list_XenaPancanExp_fn):
    with open(gene_list_XenaPancanExp_fn,'r') as gene_list_XenaPancanExp_f:
        g_list = json.load(gene_list_XenaPancanExp_f)
    return g_list

def prepare_dict_of_randGeneSets(gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp,num_of_sets_in_oneRead):

    comon_genes_XenaPancan_TCGAbrca_s = prep_common_genes_set(gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp)

    print 'Final set length: ', len(comon_genes_XenaPancan_TCGAbrca_s)

    dict_of_genel_4run={}
    for n_run in range(num_of_sets_in_oneRead):
        gene_l = list(rd.sample(comon_genes_XenaPancan_TCGAbrca_s, 39))
        dict_of_genel_4run[n_run]= gene_l
    return dict_of_genel_4run

def prep_common_genes_set(gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp):
    gene_list_XenaPancanExp_s = set(gene_list_XenaPancanExp)
    gene_list_XenaTCGAbrcaExp_s = set(gene_list_XenaTCGAbrcaExp)

    common_genes_set = gene_list_XenaPancanExp_s.intersection(gene_list_XenaTCGAbrcaExp_s)

    return common_genes_set

def prepare_TCGAGTEXallcan_df(norm_count_fn,gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp):
    common_genes_list = list(prep_common_genes_set(gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp))
    print 'len of common gene list ', str(len(common_genes_list))

    TCGAGTEX_allcan_allcol_df = pd.read_csv(norm_count_fn,header=0, index_col=0, compression='gzip',sep='\t')
    TCGAGTEX_allcan_df = TCGAGTEX_allcan_allcol_df.loc[common_genes_list,:]

    cn = list(TCGAGTEX_allcan_df.columns.values)
    rn = list(TCGAGTEX_allcan_df.index)
    print 'columnAll:', str(len(cn)), ' rowsAll: ', str(len(rn))

    return TCGAGTEX_allcan_df

def read_only_HOX_genes(hox_gene_list_fn):
    hox_g_list = []
    with open(hox_gene_list_fn,'r') as hox_gene_list_f:
        for hox_g in hox_gene_list_f:
            gene_name = hox_g.rstrip('\r\n')
            hox_g_list.append(gene_name)
    print 'no of hox_glist: ',str(len(hox_g_list)), hox_g_list
    dict_of_genel_4run={}
    dict_of_genel_4run['HOX'] = hox_g_list
    return dict_of_genel_4run


def read_file4randomgenes(norm_count_fn,tcga_gtex_exp_f,rand_genes_l,rand_selected_glist_fn):

    rand_selected_glist_f = open(rand_selected_glist_fn,'w')
    hox_exp_dict = {}
    g_list = []
    with gzip.open(norm_count_fn, 'rb') as norm_count_f:
#        ln = 0
        for gline in norm_count_f:
            line_l = gline.split('\t')
            if line_l[0] == 'sample':
                tcga_gtex_exp_f.write(gline)
                samp_list = prepare_sample_list(line_l)
            else:
                gname = line_l[0]
                if gname in rand_genes_l:
#                    print gname
                    g_list.append(gname)
                    glist_line4print = gname + '\n'
                    rand_selected_glist_f.write(glist_line4print)
                    tcga_gtex_exp_f.write(gline)

                    for ind in range(1,len(line_l)):
                        g_exp= line_l[ind].rstrip()
                        hox_exp_dict.setdefault(gname,[]).append(g_exp)

#            ln += 1
    rand_selected_glist_f.close()
#    print_dict_keys4debug(hox_exp_dict,samp_list)

    return samp_list,hox_exp_dict, g_list


def prepare_sample_list(line_l):

    samp_list =[]
    for ind in range (1, len(line_l)):
        samp_id = line_l[ind].strip()
        samp_list.append(samp_id)

    return samp_list

def print_dict_keys4debug(hox_exp_dict,samp_list):
    ind = 0
    for k in hox_exp_dict.keys():
        ind +=1
        print str(ind), k
    print '==> len sample_list:', str(len(samp_list))

    return

def prepare_TCGA_GTEX_samp_list(samp_list):

    TCGA_GTEX_samp_list = []
    for sample in samp_list:
        if (sample[:4] == 'TCGA' or sample[:4]=='GTEX'):
            if sample not in TCGA_GTEX_samp_list:
                TCGA_GTEX_samp_list.append(sample)
    print 'all samp len:', str(len(samp_list)), 'TCGA_GTEX len: ', str(len(TCGA_GTEX_samp_list))

    return TCGA_GTEX_samp_list

def prepare_category_df(TCGA_GTEX_samp_list,samp_category_dict,samp_det_category_dict):

    sample_category_dict4df ={}
    for sample in TCGA_GTEX_samp_list:

        if sample in samp_category_dict:
            samp_category = samp_category_dict[sample][0]
        else:
            samp_category = 'NoCategory'
        sample_category_dict4df.setdefault(sample,[]).append(samp_category)
        if sample in samp_det_category_dict:
            samp_det_cat = samp_det_category_dict[sample][0]
        else:
            samp_det_cat = 'NoDetCat'
        sample_category_dict4df.setdefault(sample,[]).append(samp_det_cat)

#    line4prnt,cohort = prepare_line4print(sample, samp_category, samp_det_cat, hox_genes_l, HOX_dict, samp_no)

    return sample_category_dict4df



def prepare_dfW100randSamples_of_every_cohort(TCGA_GTEX_breast_df,expression_file_source):

    tcga_allBreast_rownames = list(TCGA_GTEX_breast_df.index)
    if expression_file_source == 'TCGAonly':
        TCGA_breast_tumor_sampl, GTEX_breast_health_sampl = prep_list_of_tumorHealth_TCGAonly(tcga_allBreast_rownames)
    else:
        if expression_file_source == 'TCGA_GTEX':
            TCGA_breast_tumor_sampl, GTEX_breast_health_sampl = prep_list_of_tumorHealth_TCGA_GTEX(tcga_allBreast_rownames)
        else:
            print ' ===> Error Error invalid expression file source', expression_file_source

# At this point there are 2 lists of samples GTEX (for health sample IDs) and TCGA (for tumor sample IDs)
# from every list 100 samples are randomly selected (to compare)
#    print sample(TCGA_breast_tumor_sampl)
###### Need to remove health samples from
    TCGA_breast_tumor_100_randSampl = sample_list(TCGA_breast_tumor_sampl, 100)
    GTEX_breast_health_100_randSampl = sample_list(GTEX_breast_health_sampl, 100)
    TCGA_GTEX_breast_200randSamples_l = join2lists(TCGA_breast_tumor_100_randSampl,GTEX_breast_health_100_randSampl)
    TCGA_GTEX_breast_200randSamples_df = TCGA_GTEX_breast_df.loc[TCGA_GTEX_breast_200randSamples_l]


    return TCGA_GTEX_breast_200randSamples_df, TCGA_breast_tumor_100_randSampl

def prepare_dfWTCGA_N100randSamples_of_GTEX(TCGA_GTEX_breast_df,expression_file_source,TCGAonly_tumor_samp_l):

    tcga_allBreast_rownames = list(TCGA_GTEX_breast_df.index)
    if expression_file_source == 'TCGA_GTEX':
        TCGA_breast_tumor_sampl, GTEX_breast_health_sampl = prep_list_of_tumorHealth_TCGA_GTEX(tcga_allBreast_rownames)
    else:
        print ' ===> Error Error invalid expression file source', expression_file_source

# At this point there are 2 lists of samples GTEX (for health sample IDs) and TCGA (for tumor sample IDs)
# we would like to select the same 100 samples form TCGA only list and 100 random samples from GTEX (to compare)
#
    GTEX_breast_health_100_randSampl = sample_list(GTEX_breast_health_sampl, 100)
    TCGA_GTEX_breast_200randSamples_l = join2lists(TCGAonly_tumor_samp_l,GTEX_breast_health_100_randSampl)
    TCGA_GTEX_breast_200randSamples_df = TCGA_GTEX_breast_df.loc[TCGA_GTEX_breast_200randSamples_l]

    return TCGA_GTEX_breast_200randSamples_df


def sample_list(list2sample,num_of_samples):

    sampled_list_of_obj = []
    length_oflist = len(list2sample) - 1
    num_l = rd.sample(range(0,length_oflist),num_of_samples)
    for rnd_num in num_l:
        sampled_list_of_obj.append(list2sample[rnd_num])
    print 'length of sampled list is: ', str(len(sampled_list_of_obj)), str(length_oflist)

    return sampled_list_of_obj

def join2lists(TCGA_breast_tumor_100_randSampl,GTEX_breast_health_100_randSampl):
    oneBiglist = []
    map(oneBiglist.extend, [TCGA_breast_tumor_100_randSampl,GTEX_breast_health_100_randSampl])
    return oneBiglist


def prepare_breast_TCGA_GTEX_df4RandGenesN(rand_Glist,TCGA_GTEX_df,TCGAonly_tumor_samp_l):

    allSamp_rand_gene_df = TCGA_GTEX_df.loc[randg_list,:]
#    allSamp_rand_gene_df.columns = samp_list
#    allSamp_rand_gene_df.to_csv(debug_df_fn,sep=',',index=True)  # to debug the created df
    samp_list = list(allSamp_rand_gene_df.columns.values)
    gl = list(allSamp_rand_gene_df.index)                         # to debug
# Since the GTEX_TCGA file includes also Target samples an K samples, so need to remove them
    TCGA_GTEX_samp_list = prepare_TCGA_GTEX_samp_list(samp_list)
    print 'all samp len:', str(len(samp_list)), 'TCGA_GTEX len: ', str(len(TCGA_GTEX_samp_list)), 'num of genes in allsamp_df: ', str(len(gl))
    TCGA_GTEX_rand_gene_df = allSamp_rand_gene_df[TCGA_GTEX_samp_list]    # make sure df includes only tcga & gtex samples
    TCGA_GTEX_rand_gene_traspose_df = TCGA_GTEX_rand_gene_df.transpose()  # now column names=genes rownames=samples

#    TCGA_GTEX_rand_gene_traspose_df.to_csv(debug_dft_fn,sep=',',index=True)  # to debug the created df

    sample_category_dict4df = prepare_category_df(TCGA_GTEX_samp_list,samp_category_dict,samp_det_category_dict)

    samp_category_df = pd.DataFrame.from_dict(sample_category_dict4df,orient='index')
    samp_category_df.columns = ['category','detailed_category']
    TCGA_GTEX_rand_genWcategory_df = pd.concat([TCGA_GTEX_rand_gene_traspose_df,samp_category_df],axis=1,join_axes=[TCGA_GTEX_rand_gene_traspose_df.index])

#    TCGA_GTEX_rand_genWcategory_df.to_csv(debug_dftcat_fn,sep=',',index=True)  # to debug the created df

    print 'len of sample_category_dict4df:', str(len(sample_category_dict4df))
##pprint.pprint(sample_category_dict4df, width=1)

    breast_category_l = ['TCGA Breast Invasive Carcinoma','GTEX Breast']
    TCGA_GTEX_breast_df = TCGA_GTEX_rand_genWcategory_df.loc[TCGA_GTEX_rand_genWcategory_df['category'].isin(breast_category_l)]
#    TCGA_GTEX_breast_df.to_csv(debug_dfbreast_fn,sep=',',index=True)   # to debug the created df

    expression_file_source = 'TCGA_GTEX'
    TCGA_GTEX_breast_rand100Tn100H_df = prepare_dfWTCGA_N100randSamples_of_GTEX(TCGA_GTEX_breast_df,expression_file_source,TCGAonly_tumor_samp_l)
    TCGA_GTEX_breast_rand100Tn100H_df = TCGA_GTEX_breast_rand100Tn100H_df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')

    return TCGA_GTEX_breast_rand100Tn100H_df

def check_Wilcox4both_source(TCGAGTEX_200sampBreast_df,TCGAonly_200sampBreast_df,Wilcox_output_dir,run_no,summary_Wilcox_out_f):
    TCGAonly_col_name_l = list(TCGAonly_200sampBreast_df.columns.values)
    TCGAGTEX_col_name_l = list(TCGAGTEX_200sampBreast_df.columns.values)

    TCGAonly_row_name_l = list(TCGAonly_200sampBreast_df.index)
    TCGAGTEX_row_name_l = list(TCGAGTEX_200sampBreast_df.index)

    common_genes_set = prep_common_genes_set(TCGAonly_col_name_l,TCGAGTEX_col_name_l)
    print '\n len of common gene names:', str(len(list(common_genes_set)))

    TCGAonly_YES=0
    TCGAonly_NO=0
    TCGAGTEX_YES=0
    TCGAGTEX_NO=0
    spDiff=0
    spDiff2fold=0
    diff_gene_list =[]
    diff_gene_list2fold =[]
    for col_nm in TCGAonly_col_name_l:
        TCGAonly_col_nm_tumor_exp_l, TCGAonly_col_nm_health_exp_l = prep_list_of_tumorHealth_TCGAonly(TCGAonly_row_name_l)
        if col_nm in TCGAGTEX_col_name_l:
            TCGAGTEX_col_nm_tumor_exp_l, TCGAGTEX_col_nm_health_exp_l = prep_list_of_tumorHealth_TCGA_GTEX(TCGAGTEX_row_name_l)

            one_gene_TCGAonly_tumor = list(TCGAonly_200sampBreast_df[col_nm].loc[TCGAonly_col_nm_tumor_exp_l])
            one_gene_TCGAonly_health = list(TCGAonly_200sampBreast_df[col_nm].loc[TCGAonly_col_nm_health_exp_l])

            syg_yesTCGA, no_sygTCGA, GT2_TCGA = run_Wilcox_one_gene(col_nm,one_gene_TCGAonly_tumor,one_gene_TCGAonly_health,Wilcox_output_f,'TCGAonly',run_no)
            TCGAonly_NO = TCGAonly_NO + no_sygTCGA
            TCGAonly_YES = TCGAonly_YES + syg_yesTCGA

            one_gene_TCGAGTEX_tumor = list(TCGAGTEX_200sampBreast_df[col_nm].loc[TCGAGTEX_col_nm_tumor_exp_l])
            one_gene_TCGAGTEX_health = list(TCGAGTEX_200sampBreast_df[col_nm].loc[TCGAGTEX_col_nm_health_exp_l])

            syg_yesTCGAGTEX, no_sygTCGAGTEX, GT2_TCGAGTEX = run_Wilcox_one_gene(col_nm,one_gene_TCGAGTEX_tumor,one_gene_TCGAGTEX_health,Wilcox_output_f,'TCGAGTEX',run_no)
            TCGAGTEX_NO = TCGAGTEX_NO + no_sygTCGAGTEX
            TCGAGTEX_YES = TCGAGTEX_YES + syg_yesTCGAGTEX

            if (syg_yesTCGA==1 and GT2_TCGA==1):
                if (syg_yesTCGAGTEX==1 and GT2_TCGA==1):
                    continue
                else:
                    spDiff2fold += 1
                    diff_gene_list2fold.append(col_nm)
            else:
                if (syg_yesTCGAGTEX==1 and GT2_TCGAGTEX==1):
                    spDiff2fold += 1
                    diff_gene_list2fold.append(col_nm)

            if syg_yesTCGA != syg_yesTCGAGTEX:
                spDiff +=1
                diff_gene_list.append(col_nm)

    dif = TCGAonly_YES - TCGAGTEX_YES
    gene_list_str = " ".join(str(x) for x in TCGAonly_col_name_l)
    difg_list_str = " ".join(str(x) for x in diff_gene_list)
    difg_list2fold_str = " ".join(str(x) for x in diff_gene_list2fold)
    sum_line = (str(run_no) + ',' + str(TCGAonly_YES) + ',' + str(TCGAonly_NO) + ',' + str(TCGAGTEX_YES) + ',' +
                str(TCGAGTEX_NO)  + ',' + str(dif) + ',' + str(spDiff) + ',' + str(spDiff2fold) + ',' +
                difg_list_str + ','+ difg_list2fold_str + ','+ gene_list_str + '\n')
    summary_Wilcox_out_f.write(sum_line)


    return

def run_Wilcox_one_gene(gene_nm,one_gene_TCGAonly_tumor,one_gene_TCGAonly_health,Wilcox_output_f,exp_source,run_no):

    one_gene_TCGAonly_tumor_f= np.array(one_gene_TCGAonly_tumor).astype(np.float)
    one_gene_TCGAonly_health_f = np.array(one_gene_TCGAonly_health).astype(np.float)
    one_gene_TCGAonly_Wilcox = sc.ranksums(one_gene_TCGAonly_tumor_f,one_gene_TCGAonly_health_f)

#    print '\n',exp_source,' Wilcox4 gene: ', gene_nm, one_gene_TCGAonly_Wilcox
    line4print = (str(run_no) + ',' + gene_nm + ',' + exp_source + ',' + str(np.mean(one_gene_TCGAonly_health_f)) + ',' +
                  str(np.mean(one_gene_TCGAonly_tumor_f)) + ',' +
                  str(one_gene_TCGAonly_Wilcox.statistic) + ',' +
                 str(one_gene_TCGAonly_Wilcox.pvalue))
    adjust_pval = 0.05/39 #bonferroni adjusted p_value

   # if one_gene_TCGAonly_Wilcox.pvalue < 0.05:
    if one_gene_TCGAonly_Wilcox.pvalue < adjust_pval:
        line4print = line4print + ',YES'
        syg_yes = 1
        no_syg = 0
    else:
        line4print = line4print + ',NO'
        syg_yes = 0
        no_syg = 1

    if np.mean(one_gene_TCGAonly_health_f) != 0:
        exp_change = float(np.mean(one_gene_TCGAonly_tumor_f)/np.mean(one_gene_TCGAonly_health_f))
        if (exp_change>2 or exp_change<0.5):
            line4print = line4print + ',' + str(exp_change) + ',' + 'GT-2\n'
            is_exp_gt2 = 1
        else:
            line4print = line4print + ',' + str(exp_change) + ',' + 'NO\n'
            is_exp_gt2 = 0
    else:
        if np.mean(one_gene_TCGAonly_tumor_f) != 0:
            line4print = line4print + ',' + str(np.mean(one_gene_TCGAonly_tumor_f)) + ',' + 'GT-0\n'
            is_exp_gt2 = 0
        else:
            line4print = line4print + ',' + str(0) + ',' + 'NO\n'
            is_exp_gt2 = 0

    Wilcox_output_f.writelines(line4print)

    return syg_yes, no_syg, is_exp_gt2

source_dir = "/home/HOX/Xena"
category_tab_fn = os.path.join(source_dir,"TCGA_GTEX_category.txt","TCGA_GTEX_category.txt")    # downloaded from Xena
det_cat_tab_fn = os.path.join(source_dir,"TCGA_GTEX_category.txt","TcgaTargetGTEX_selected_phenotype.txt")  # downloaded from Xena

norm_cFromZena_fn = "TcgaTargetGtex_RSEM_Hugo_norm_count.gz"    # downloaded from Xena
TCGA_BRCA_exp_fromXena_fn = os.path.join(source_dir,"TCGA_BRCAexpression","HiSeqV2.gz")   # downloaded from Xena

gene_list_XenaPancanExp_fn = os.path.join(source_dir,"testRandom","gene_list_XenaPancanExp.json")
gene_list_XenaTCGAbrcaExp_fn = os.path.join(source_dir,"testRandom","gene_list_XenaTCGAbrcaExp.json")

hox_gene_list_fn = "/home/HOX/Glist/glist39HOX.csv"

### Output files

tcga_gtex_exp_fn = os.path.join(source_dir,"testRandom","TCGA_GTEX_exp_norm_count.csv")
tcga_brca_xena_exp_fn = os.path.join(source_dir,"testRandom","TCGA_brca__xena_exp_norm.csv")
rand_selected_glist_dir = os.path.join(source_dir,"testRandom/Glist","")
final_rangenes_df_dir = os.path.join(source_dir,"testRandom/randgenesRuns","")          # for debug
debug_TCGAonly_dfbreast_fn = os.path.join(source_dir,"testRandom","TCGAonly_randgene1_dfbreast.csv") # for debug
debug_TCGAGTEXsamp_dfbreast_fn = os.path.join(source_dir,"testRandom","TCGAGTEX_randgene1_dfbreast.csv")

#Wilcox_output_fn = "/home/HOX/Xena/testRandom/ttest/20200122/Wilcox_results_randgene1_breast.csv"
Wilcox_output_dir = os.path.join(source_dir,"testRandom/ttest/20200123","")

###
###   Begining of the program
###
###
now = datetime.now()
print 'start time: ', now.isoformat()
samp_category_dict =  prep_sample_category_dict(category_tab_fn)
samp_det_category_dict =  prep_sample_category_dict(det_cat_tab_fn)

norm_count_fn = os.path.join(source_dir,norm_cFromZena_fn)

num_of_sets_in_oneRead = 1000
# num_of_sets_in_oneRead = 10   # for debug
runoutput_fn = 'allRunSum_Wilcox_results.csv'

print 'num of renadom sets in this run is:', str(num_of_sets_in_oneRead)
gene_list_XenaPancanExp = load_json2list(gene_list_XenaPancanExp_fn)

gene_list_XenaTCGAbrcaExp = load_json2list(gene_list_XenaTCGAbrcaExp_fn)

dict_of_genel_4run= prepare_dict_of_randGeneSets(gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp,num_of_sets_in_oneRead)

#dict_of_genel_4run = read_only_HOX_genes(hox_gene_list_fn)    # For HOX only

##pprint.pprint (dict_of_genel_4run,width=1)

### prepare the df that includes the expression data of BRCA samples from TCGA only
TCGAonly_BRCA_df = pd.read_csv(TCGA_BRCA_exp_fromXena_fn,header=0, index_col=0, compression='gzip',sep='\t')

samples_notIn_TCGA_GTEX =['TCGA-AR-A0U1-01','TCGA-BH-A0AY-01','TCGA-C8-A1HF-01','TCGA-E2-A108-01','TCGA-E2-A1IP-01',
                          'TCGA-BH-A0B2-11']
TCGAonly_BRCA_df.drop(samples_notIn_TCGA_GTEX,inplace=True,axis=1)

#### prepare the df that includes the expression data of TCGAGTEX 4 all cancers
TCGAGTEX_allcan_df = prepare_TCGAGTEXallcan_df(norm_count_fn,gene_list_XenaPancanExp,gene_list_XenaTCGAbrcaExp)



# prepare the summary Wilcox result file
summary_Wilcox_out_fn = os.path.join(Wilcox_output_dir, runoutput_fn)
summary_Wilcox_out_f = open(summary_Wilcox_out_fn,'w')
title_sum = 'run_no,TCGA_only-YES,TCGA_only-No,TCGAGTEX-Yes,TCGAGTEX-No,Diff,SpDiff,SpDiff2fold,DiffL,Diff2foldL,GeneList\n'
summary_Wilcox_out_f.write(title_sum)

Wilcox_output_fn = os.path.join(Wilcox_output_dir, 'runNoAll' + '_Wilcox_results_randgenes_breast.csv')
Wilcox_output_f = open(Wilcox_output_fn,'w')
title = 'runno,gene,Exp_source,Health_mean,Tumor_mean,statistics,pval,is_syg,exp_change,GT2\n'
Wilcox_output_f.write(title)


for r_list in dict_of_genel_4run.keys():
    randg_list = dict_of_genel_4run[r_list]
    print 'runno: ', str(r_list), randg_list

# # prepare 2 gene lists:
# # 1. list of genes that appear in the XenaPancan expression file
# # 2. list of the genes that appear in the Xena TCGAbrca expression file
# #
# gene_list_XenaPancanExp = prepare_gene_list(norm_count_fn)
# write_l2_json(gene_list_XenaPancanExp, gene_list_XenaPancanExp_fn)
#
# gene_list_XenaTCGAbrcaExp = prepare_gene_list(TCGA_BRCA_exp_fromXena_fn)
# write_l2_json(gene_list_XenaTCGAbrcaExp, gene_list_XenaTCGAbrcaExp_fn)

#################################################################################################################
# Randomly select genes from the list that appears in the the BRCA expression file that was downloaded from Xena
#
    list_num = r_list
    run_no = r_list
    rand_selected_glist_fn = rand_selected_glist_dir + 'glist' + str(list_num) + '.txt'

    ###############################################################################################################
    ## after the command below:
    ## samp_list - is a list of all the samples (actually the titles) for which expression exist in the file)
    ## hox_ecp_dict - is a dictionary where key is the selected genename and the value is a list of expression of all
    ##                the samples that appear in samp_list
    # tcga_brca_xena_exp_f = open(tcga_brca_xena_exp_fn,'w')    # for debug purposes only
    # samp_list,hox_exp_dict,rand_Glist = read_file4randomgenes(TCGA_BRCA_exp_fromXena_fn,tcga_brca_xena_exp_f,ln_list,rand_selected_glist_fn)
    # tcga_brca_xena_exp_f.close()
    #
    # # build a df from the dictionary that include the genes and their lists of expressions
    # #
    # rand_gene_df = pd.DataFrame.from_dict(hox_exp_dict,orient='index')
    # rand_gene_df.columns = samp_list
    # rand_gene_df.to_csv(debug_df_fn,sep=',',index=True)  # to debug the created df
    #
    # print 'all TCGA BRCA tumor health samp len:', str(len(samp_list)),'samp_list len: ', str(len(samp_list))
    rand_gene_df = TCGAonly_BRCA_df.loc[randg_list,:]
    cn = list(rand_gene_df.columns.values)
    rn = list(rand_gene_df.index)
    print 'column:', str(len(cn)), ' rows: ', str(len(rn))
    print 'rn', rn

    TCGA_brca_rand_gene_traspose_df = rand_gene_df.transpose()  # now column names=genes rownames=samples
#    TCGA_brca_rand_gene_traspose_df.to_csv(debug_brca_dft_fn,sep=',',index=True)  # to debug the created df


    expression_file_source = 'TCGAonly'
    TCGA_breast_rand100Tn100H_df,TCGAonly_100rand_tumorSamp_l = prepare_dfW100randSamples_of_every_cohort(TCGA_brca_rand_gene_traspose_df,expression_file_source)
    TCGA_breast_rand100Tn100H_df = TCGA_breast_rand100Tn100H_df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')

    debug_TCGAonly_dfbreast_fn = final_rangenes_df_dir + 'run' + str(run_no) + '_TCGAonly_randg_breastDf.csv'
    TCGA_breast_rand100Tn100H_df.to_csv(debug_TCGAonly_dfbreast_fn,sep=',',index=True)   # to debug the created df

    # at this point we have a df based on TCGAonly expression that includes 100 health samples and 100 tumor samples
    # (randomly slected) on 39 randomly selected genes.
    # Now we want to prepare a df based on TCGA_GTEX that will include same genes, w same tumor samples but different health
    # samples

    TCGA_GTEX_breast_rand100Tn100H_df = prepare_breast_TCGA_GTEX_df4RandGenesN(randg_list,TCGAGTEX_allcan_df,TCGAonly_100rand_tumorSamp_l)
    debug_TCGAGTEXsamp_dfbreast_fn = final_rangenes_df_dir + 'run' + str(run_no) + '_TCGAGTEX_randg_breastDf.csv'
    TCGA_GTEX_breast_rand100Tn100H_df.to_csv(debug_TCGAGTEXsamp_dfbreast_fn,sep=',',index=True)   # to debug the created df

    check_Wilcox4both_source(TCGA_GTEX_breast_rand100Tn100H_df,TCGA_breast_rand100Tn100H_df,Wilcox_output_dir,run_no,summary_Wilcox_out_f)

Wilcox_output_f.close()
summary_Wilcox_out_f.close()
now = datetime.now()
print 'end time: ', now.isoformat()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

