# created on 26/12 in order to run the ttest only once on all samples for every cancer type
# change on 01/01/2020 - replace ttest by wilcoxon rank test
from hoxcancer_functions import *
import os.path
import random as rd
import pandas as pd
import json
from scipy import stats as sc
import numpy as np
import glob
import pprint
import collections
#
#  This program checks for all HOX genes in all cancer types whether they are differentially expressed
#
#

def  get_cancerType_from_fn(base_fn):
    right_fn = base_fn[10:]
    cancer_type_l = right_fn.split('.')
    cancer_type = cancer_type_l[0]
    return cancer_type


def choose_random_100h100t(tumor_sample_list_cleanFromH,GTEX_sample_list):

    tumor_s = set(tumor_sample_list_cleanFromH)
    health_s = set(GTEX_sample_list)
    tumor100_samples_l = list(rd.sample(tumor_s,100))
    health100_samples_l = list(rd.sample(health_s,100))

    return health100_samples_l, tumor100_samples_l

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


def print_dict_keys4debug(hox_exp_dict,samp_list):
    ind = 0
    for k in hox_exp_dict.keys():
        ind +=1
        print str(ind), k
    print '==> len sample_list:', str(len(samp_list))

    return


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

def prepare_df4_everyCan_Nwilcox(TCGA_GTEX_allcan_df,cancer_type_samplel_dict,canType_name):

    allcantype_rownames_l = cancer_type_samplel_dict[canType_name]
    TCGA_GTEX_brain_200randSamples_df = TCGA_GTEX_allcan_df.loc[allcantype_rownames_l]
    TCGA_GTEX_brain_200randSamples_df = TCGA_GTEX_brain_200randSamples_df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')

    sampOf200_l = list(TCGA_GTEX_brain_200randSamples_df.index)
    TCGAGTEX_col_name_l=list(TCGA_GTEX_brain_200randSamples_df.columns.values)
#    print '==> cantype:', canType_name,'len of columns:', len(TCGAGTEX_col_name_l)

    return TCGA_GTEX_brain_200randSamples_df, sampOf200_l

def print_cancer_df_cleanFromH(cancer_type,tumor_sample_list_cleanFromH,GTEX_sample_list,cancerType_39HOX_df,outdir_woTCGAhealth):

    canType_FullSampl_wo_TCGAh_l=join2lists(tumor_sample_list_cleanFromH,GTEX_sample_list)
    cancerType_cleanH_fullsamples_df = cancerType_39HOX_df.loc[canType_FullSampl_wo_TCGAh_l]
    cantype_fn=cancer_type + '_woTCGAh.csv'
    out_fn = os.path.join(outdir_woTCGAhealth,cantype_fn)
    cancerType_cleanH_fullsamples_df.to_csv(out_fn,sep=',',index=True)
    print cancer_type," tumor: ", str(len(tumor_sample_list_cleanFromH))," GTEX:", str(len(GTEX_sample_list))," all", str(len(canType_FullSampl_wo_TCGAh_l))

    return

def check_wilcox4canTypeHOX(TCGA_GTEX_cantype_200randSamples_df,wilcox_output_dir,run_no,canType_sum_dict,canType_name,smpl_num):

    TCGAGTEX_col_name_l = list(TCGA_GTEX_cantype_200randSamples_df.columns.values)
    TCGAGTEX_row_name_l = list(TCGA_GTEX_cantype_200randSamples_df.index)

    TCGAGTEX_YES=0
    both_pvalMinus10Ngt2=0
    min10=0
    exp_change2zero = 0

    for col_nm in TCGAGTEX_col_name_l:
        TCGAGTEX_col_nm_tumor_exp_l, TCGAGTEX_col_nm_health_exp_l = prep_list_of_tumorHealth_TCGA_GTEX(TCGAGTEX_row_name_l)
        one_gene_TCGAGTEX_tumor = list(TCGA_GTEX_cantype_200randSamples_df[col_nm].loc[TCGAGTEX_col_nm_tumor_exp_l])
        one_gene_TCGAGTEX_health = list(TCGA_GTEX_cantype_200randSamples_df[col_nm].loc[TCGAGTEX_col_nm_health_exp_l])
        syg_yes, min10syg_yes, both_cond, change2zero = run_wilcox_one_geneHOX(col_nm,one_gene_TCGAGTEX_tumor,one_gene_TCGAGTEX_health,wilcox_output_f,canType_name,smpl_num)
 #       TCGAGTEX_YES = TCGAGTEX_YES + syg_yes
 #       min10 = min10 + min10syg_yes
        both_pvalMinus10Ngt2 =  both_pvalMinus10Ngt2 + both_cond
        exp_change2zero = exp_change2zero + change2zero
    canType_sum_dict[canType_name] = both_pvalMinus10Ngt2


     #       one_gene_TCGAGTEX_ttest = sc.ttest_ind(one_gene_TCGAGTEX_tumor,one_gene_TCGAGTEX_health,equal_var=True)
     #       print '\nTCGAGTEX ttest4 gene: ', col_nm, one_gene_TCGAGTEX_ttest
    # print 'change2zero: ', str(exp_change2zero), canType_name

    return exp_change2zero

def print_summary_line(summary_wilcox_out_f, TCGAGTEX_col_name_l,run_no,cantype_dict):
    gene_list_str = " ".join(str(x) for x in TCGAGTEX_col_name_l)
    sum_line = (str(run_no) + ',' + str(cantype_dict['Adrenal_Gland']) + ',' + str(cantype_dict['brain']) + ',' +
                                    str(cantype_dict['breast']) + ',' + str(cantype_dict['Colon']) + ',' +
                                    str(cantype_dict['Esophagus']) + ',' + str(cantype_dict['Leukemia']) + ',' +
                                    str(cantype_dict['Liver']) + ',' + str(cantype_dict['Lung']) + ',' +
                                    str(cantype_dict['Pancreas']) + ',' + str(cantype_dict['Prostate']) + ',' +
                                    str(cantype_dict['Stomach']) + ',' + str(cantype_dict['Thyroid']) + ',' +
                                    str(cantype_dict['change2zero']))
    totalsyg = 0
    for canname in cantype_dict.keys():
        totalsyg = totalsyg + cantype_dict[canname]
    sum_line = sum_line +',' + str(totalsyg) + ',' + gene_list_str + '\n'   # for total
    summary_wilcox_out_f.write(sum_line)
    return

def sample_list(list2sample,num_of_samples):

    sampled_list_of_obj = []
    length_oflist = len(list2sample) - 1
    num_l = rd.sample(range(0,length_oflist),num_of_samples)
    for rnd_num in num_l:
        sampled_list_of_obj.append(list2sample[rnd_num])
    print 'length of sampled list is: ', str(len(sampled_list_of_obj)), str(length_oflist)

    return sampled_list_of_obj


def run_wilcox_one_geneHOX(gene_nm,one_gene_TCGAonly_tumor,one_gene_TCGAonly_health,wilcox_output_f,canType_name,smpl_num):

    one_gene_TCGAGTEX_tumor_f= np.array(one_gene_TCGAonly_tumor).astype(np.float)
    one_gene_TCGAGTEX_health_f = np.array(one_gene_TCGAonly_health).astype(np.float)
    one_gene_TCGAGTEX_wilcox = sc.ranksums(one_gene_TCGAGTEX_tumor_f,one_gene_TCGAGTEX_health_f )
 #   one_gene_TCGAGTEX_ttest = sc.ttest_ind(one_gene_TCGAGTEX_tumor_f,one_gene_TCGAGTEX_health_f,equal_var=False )
#    print '\n',exp_source,' ttest4 gene: ', gene_nm, one_gene_TCGAonly_ttest
    line4print = (str(smpl_num) + ',' + canType_name + ',' + gene_nm + ',' + str(np.mean(one_gene_TCGAGTEX_health_f)) + ',' +
                  str(np.mean(one_gene_TCGAGTEX_tumor_f)) + ',' + str(one_gene_TCGAGTEX_wilcox.statistic) + ',' +
                  str(one_gene_TCGAGTEX_wilcox.pvalue) + ',' + str(np.std(one_gene_TCGAGTEX_health_f)) + ',' +
                  str(np.std(one_gene_TCGAGTEX_tumor_f)))
    change2zero = 0
    adjust_pval = 0.05/39 # bonferroni adjusted p_value
    if one_gene_TCGAGTEX_wilcox.pvalue < adjust_pval:
#    if one_gene_TCGAGTEX_wilcox.pvalue < 0.05:
        line4print = line4print + ',YES'
        syg_yes = 1
    else:
        line4print = line4print + ',NO'
        syg_yes = 0
    minus10 = -10
    if one_gene_TCGAGTEX_wilcox.pvalue < 10**minus10:
        line4print = line4print + ',YES'
        min10syg_yes = 1
    else:
        line4print = line4print + ',NO'
        min10syg_yes = 0
    if np.mean(one_gene_TCGAGTEX_tumor_f) > np.mean(one_gene_TCGAGTEX_health_f):
        line4print = line4print + ',UP'
    else:
        if np.mean(one_gene_TCGAGTEX_tumor_f) < np.mean(one_gene_TCGAGTEX_health_f):
            line4print = line4print + ',DOWN'
        else:
            line4print = line4print + ',NOchange'

    if np.mean(one_gene_TCGAGTEX_health_f) != 0:
        exp_change = float(np.mean(one_gene_TCGAGTEX_tumor_f)/np.mean(one_gene_TCGAGTEX_health_f))
        if (exp_change>2 or exp_change<0.5):
            line4print = line4print + ',' + str(exp_change) + ',' + 'GT-2'
            is_exp_gt2 = 1
            if (exp_change>4 or exp_change<0.25):
                line4print = line4print + ',' + 'GT-4\n'
            else:
                line4print = line4print + ',' + 'NO\n'
        else:
            line4print = line4print + ',' + str(exp_change) + ',' + 'NO,NO\n'
            is_exp_gt2 = 0
    else:
        if np.mean(one_gene_TCGAGTEX_tumor_f) != 0:
            line4print = line4print + ',' + str(np.mean(one_gene_TCGAGTEX_tumor_f)) + ',' + 'GT-0\n'
            is_exp_gt2 = 0
            if np.mean(one_gene_TCGAGTEX_tumor_f) > 2:  # if it is significant change
                change2zero = 1
            else:
                change2zero=0
        else:
            line4print = line4print + ',' + str(0) + ',' + 'NO\n'
            is_exp_gt2 = 0


    if (min10syg_yes == 1 & is_exp_gt2 == 1):
        both_conditions = 1
    else:
        both_conditions = 0

    wilcox_output_f.writelines(line4print)

    return syg_yes, min10syg_yes, both_conditions, change2zero

def prepare_cantype_dictNwilcoxoutFileHOX(run_no,wilcox_output_dir):
    cantype_dict = {}
    cantype_dict['Adrenal_Gland'] = 0
    cantype_dict['Adrenal_Gland_PCPG'] = 0
    cantype_dict['brain'] = 0
    cantype_dict['brain_LGG'] = 0
    cantype_dict['brain_GBM'] = 0
    cantype_dict['breast'] = 0
    cantype_dict['Colon'] = 0
    cantype_dict['Esophagus'] = 0
    cantype_dict['Leukemia'] = 0
    cantype_dict['Liver'] = 0
    cantype_dict['Lung'] = 0
    cantype_dict['Lung_LUAD'] = 0
    cantype_dict['Lung_LUSC'] = 0
    cantype_dict['Pancreas'] = 0
    cantype_dict['Prostate'] = 0
    cantype_dict['Stomach'] = 0
    cantype_dict['Thyroid']= 0
    cantype_dict['change2zero'] = 0

    wilcox_fn = 'runNo' + str(run_no) + '_wilcox_results_randgenes.csv'
    wilcox_output_fn = os.path.join(wilcox_output_dir , wilcox_fn)
    wilcox_output_f = open(wilcox_output_fn,'w')
    title = 'sample_no,cancer_type,gene,mean_GTEX,mean_tcga,wilcox_statistics,pval,std_GTEX,std_TCGA,is_syg,is10min10syg,up_down,exp_diff,isGT2,isGT4\n'
    wilcox_output_f.write(title)

    return cantype_dict,wilcox_output_f

category_tab_fn = "/home/HOX/Xena/TCGA_GTEX_category.txt/TCGA_GTEX_category.txt"
det_cat_tab_fn = "/home/HOX/Xena/TCGA_GTEX_category.txt/TcgaTargetGTEX_selected_phenotype.txt"

source_dir = "/home/HOX/Xena/TwlvSamplists"
#source_dir = "/home/HOX/Xena/testRandom4allCancers/BrainSamp"  # 23/01

TCGA_BRCA_exp_fromXena_fn = "/home/HOX/Xena/TCGA_BRCAexpression/HiSeqV2.gz"  ######################

gene_list_XenaPancanExp_fn = "/home/HOX/Xena/testRandom/gene_list_XenaPancanExp.json"
gene_list_XenaTCGAbrcaExp_fn = "/home/HOX/Xena/testRandom/gene_list_XenaTCGAbrcaExp.json"

hox_gene_list_fn = "/home/HOX/Glist/glist39HOX.csv"
cancer_type_smpldict_fn = "/home/HOX/Xena/testRandom4allCancers/Samplelist200/samplist_byCancerType.json"

### Output files
#2301output_dir = "/home/HOX/Xena/testRandom4allCancers/"
outdir_woTCGAhealth = "/home/HOX/Xena/TwlvSamplistsWoTCGAhealth"

#tcga_gtex_exp_fn = "/home/HOX/Xena/testRandom/TCGA_GTEX_exp_norm_count.csv"
#tcga_brca_xena_exp_fn = "/home/HOX/Xena/testRandom/TCGA_brca__xena_exp_norm.csv"
#tcga_gtex_randgenes_exp_fn = "/home/HOX/Xena/testRandom4allCancers/TCGA_gtex__randgenes_exp_norm.csv"
#rand_selected_glist_dir = "/home/HOX/Xena/testRandom4allCancers/Glist/"            #v

wilcox_output_dir = "/home/HOX/Xena/testRandom4allCancers/wilcox240120"  # 01/01/2020
###
###   Begining of the program
###
###


dict_of_genel_4run = read_only_HOX_genes(hox_gene_list_fn)  # to prepare a list of hox genes

# prepare the summary wilcox result file
# summary_wilcox_out_fn = wilcox_output_dir + 'allRunSum_wilcox_resultsHOX.csv'
# summary_wilcox_out_fn = os.path.join(wilcox_output_dir , 'allRunSum_wilcox_resultsHOX.csv')   # 23/01
# summary_wilcox_out_f = open(summary_wilcox_out_fn,'w')
# title_sum = ('run_no,Adrenal_Gland,Adrenal_Gland_PCPG,brain,brain_LGG,brain_GBM,Breast,Colon,Esophagus,Leukemia,' +
#              'Liver,Lung,Lung_LUAD,Lung_LUSC,Pancreas,Prostate,Stomach,' +
#             'Thyroid,expChange2zero,Total4run,GeneList\n')
# summary_wilcox_out_f.write(title_sum)

for r_list in dict_of_genel_4run.keys():
    HOX_list = dict_of_genel_4run[r_list]
    print 'run no:', str(r_list), HOX_list

canType_sum_dict, wilcox_output_f = prepare_cantype_dictNwilcoxoutFileHOX('HOX1',wilcox_output_dir)
exp_change2zero = 0
for cancer_fn in glob.glob(os.path.join(source_dir, '*.csv')):
    base_fn = os.path.basename(cancer_fn)
    cancer_type = get_cancerType_from_fn(base_fn)
    cancerType_df =  pd.read_csv(cancer_fn,header=0, index_col=0)
    col2drop_l =['UQ-Samp_id','Order',' Cohort','Category','Detailed Category']
    cn = list(cancerType_df.columns.values)
    cancerType_df.drop(col2drop_l,axis=1,inplace=True)


    cancerType_39HOX_df = cancerType_df.loc[:,HOX_list]           # prepares a df that includes only hox genes
    n_col = len(list(cancerType_39HOX_df.columns.values))
    n_of_lines = len(list(cancerType_39HOX_df.index))
    print 'ct: ', cancer_type, 'No of colums: ', str(n_col), 'No of samples: ', str(n_of_lines)
    cancerType_row_names = list(cancerType_39HOX_df.index)
    TCGA_sample_list, GTEX_sample_list = prep_list_of_tumorHealth_TCGA_GTEX(cancerType_row_names)
    tumor_sample_list_cleanFromH = remove_health_samplesFrom_list(TCGA_sample_list)
    print cancer_type,str(len(TCGA_sample_list)),str(len(tumor_sample_list_cleanFromH)),str(len(GTEX_sample_list))
# The line below writes the files of every specific cancer type with our the health samples to a file:
 #   print_cancer_df_cleanFromH(cancer_type,tumor_sample_list_cleanFromH,GTEX_sample_list,cancerType_39HOX_df,outdir_woTCGAhealth)
    hNt_AllSamples_l = join2lists(GTEX_sample_list,tumor_sample_list_cleanFromH)
    cancer_typeAllSamples_df = cancerType_39HOX_df.loc[hNt_AllSamples_l]
    run_no = 1
    numOf_expChange2zero = check_wilcox4canTypeHOX(cancer_typeAllSamples_df,wilcox_output_dir,run_no,canType_sum_dict,cancer_type,run_no)

    exp_change2zero = exp_change2zero + numOf_expChange2zero
    canType_sum_dict['change2zero'] = exp_change2zero
wilcox_output_f.close()
hox_gl = list(cancer_typeAllSamples_df.columns.values)
