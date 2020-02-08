import os.path

#norm_count_fn = os.path.join(source_dir,norm_cFromZena_fn)
def prep_TCGA_GTEX_category_dict(category_fn):
    tcga_gtext_category_dict = {}
    fn_tcga_gtext_dict = {}
    gtexCat_fn_dict = {}
    first_line = True
    with open (category_fn, 'r') as category_f:
        for tcga_categoryl in category_f:
            if first_line:
                first_line = False
                continue
            categoryl_spl = tcga_categoryl.split(',')
            fn = categoryl_spl[2]
            det_category = categoryl_spl[0]
            cohort = categoryl_spl[1]
            tcga_gtext_category_dict.setdefault(det_category,[]).append(cohort)
            tcga_gtext_category_dict.setdefault(det_category ,[]).append(fn)

            fn_tcga_gtext_dict.setdefault(fn,[]).append(det_category )                
               
    return tcga_gtext_category_dict, fn_tcga_gtext_dict

def prep_TCGA_GTEX_category_dict_old(category_fn):
    tcga_gtext_category_dict = {}
    fn_tcga_gtext_dict = {}
    gtexCat_fn_dict = {}
    first_line = True
    with open (category_fn, 'r') as category_f:
        for tcga_categoryl in category_f:
            if first_line:
                first_line = False
                continue
            categoryl_spl = tcga_categoryl.split(',')
            fn = categoryl_spl[2].rstrip()
            tcga_cat = categoryl_spl[0]
            gtex_cat = categoryl_spl[1]
            tcga_gtext_category_dict.setdefault(tcga_cat,[]).append(gtex_cat)
            tcga_gtext_category_dict.setdefault(tcga_cat,[]).append(fn)
            if fn in fn_tcga_gtext_dict.keys():
                if tcga_cat in fn_tcga_gtext_dict[fn]:
                    pass
                else:
                    fn_tcga_gtext_dict.setdefault(fn,[]).append(tcga_cat)
            else:
                fn_tcga_gtext_dict.setdefault(fn,[]).append(tcga_cat)
                
            if gtex_cat in gtexCat_fn_dict.keys():
                if fn in gtexCat_fn_dict[gtex_cat]:
                    pass
                else:
                    gtexCat_fn_dict[gtex_cat] = fn
            else:
                gtexCat_fn_dict[gtex_cat] = fn
                
    return tcga_gtext_category_dict, fn_tcga_gtext_dict, gtexCat_fn_dict

#def read_smaple_exp_fileN_split(samplexp_fn, tcga_gtext_category_dict, fn_tcga_gtext_dict, gtexCat_fn_dict):
def read_smaple_exp_fileN_split(samplexp_fn, tcga_gtext_category_dict, fn_tcga_gtext_dict):
    first_line = True
    cohort = 'TCGA'
#    cancer_file = 
    with open(samplexp_fn,'r') as samplexp_f:
        prev_category = ' '
        for samp_line in samplexp_f:
            if first_line :
                title_line = samp_line
                first_line = False
                fn = 'TCGA_GTEX_' + 'Gtex-NoTCGA.csv'
                cancert_fn = os.path.join(source_dirn,fn) 
                cancer_file = open(cancert_fn,'w')
                cancer_file.write(title_line)
                continue
            samp_line_splt = samp_line.split(',')
            cohort = samp_line_splt[3]
            category = samp_line_splt[4]
            detailed_category = samp_line_splt[5]
            if (cohort == 'TARGET' or cohort == 'K'):
                continue
            if detailed_category != prev_category:             
#                if category != prev_category:
                cancer_file.close()
                print "category:", category, detailed_category
                # if category == 'NoCategory':
                #     print 'noCategory: ', samp_line
                #     continue
                if category != 'NoCategory':
                    fn = 'TCGA_GTEX_' + tcga_gtext_category_dict[detailed_category][1] + '.csv'
                else:
                    if detailed_category in tcga_gtext_category_dict.keys():
                        fn = 'TCGA_GTEX_' + tcga_gtext_category_dict[detailed_category][1] + '.csv'
                    else:
                        fn = 'TCGA_GTEX_NoCategory' + '.csv'
                cancert_fn = os.path.join(source_dirn,fn)
                file_exists = os.path.isfile(cancert_fn)
                if file_exists:
                    cancer_file = open(cancert_fn,'a')
                else:
                    cancer_file = open(cancert_fn,'w') 
                    cancer_file.write(title_line)
            cancer_file.write(samp_line)
            prev_category = detailed_category
    cancer_file.close()
    return


#5source_dirn = "/home/HOX/Xena/SamplistFPKM"
source_dirn = "/home/HOX/Xena/Samplists"
#1 source_dirn = "/home/HOX/Xena/SamplistZFRand"
#2 source_dirn = "/home/HOX/Xena/SamplistRBRand"
#6 source_dirn = "/home/HOX/Xena/SamplistRBslRand"
category_fn = "/home/HOX/Xena/TCGA_GTEX_category.txt/TCGA_GTEX_det_categoryl.csv"
samplexp_fn = "/home/HOX/Xena/Samplists/All/TCGA_GTEX_exp_norm_count-trasposed-Norder.csv"
#5samplexp_fn = "/home/HOX/Xena/SamplistFPKM/All/TCGA_GTEX_exp_rsem_gene_fpkm-trasposed-Nordered.csv"

#1 samplexp_fn = "/home/HOX/Xena/SamplistZFRand/All/TCGA_GTEX_exp_HUGO_norm_count-RANDZFgenes-trasposed-Ordered.csv"
#2 samplexp_fn = "/home/HOX/Xena/SamplistRBRand/All/TCGA_GTEX_exp_HUGO_norm_count-RANDRBgenes-trasposed-Ordered.csv"
#6 samplexp_fn = "/home/HOX/Xena/SamplistRBslRand/All/TCGA_GTEX_exp_HUGO_norm_count-RANDRBslgenes-trasposed-Ordered.csv"


#tcga_gtext_category_dict, fn_tcga_gtext_dict, gtexCat_fn_dict = prep_TCGA_GTEX_category_dict(category_fn)
tcga_gtext_category_dict, fn_tcga_gtext_dict = prep_TCGA_GTEX_category_dict(category_fn)

for f in fn_tcga_gtext_dict.keys():
    print '\n',f
    for i in range (len(fn_tcga_gtext_dict[f])):
        print fn_tcga_gtext_dict[f][i]

for f in tcga_gtext_category_dict.keys():
#    print '\n',f
#    for i in range (len(tcga_gtext_category_dict[f])):
    print '\n',f, tcga_gtext_category_dict[f]       

#read_smaple_exp_fileN_split(samplexp_fn, tcga_gtext_category_dict, fn_tcga_gtext_dict, gtexCat_fn_dict)
read_smaple_exp_fileN_split(samplexp_fn, tcga_gtext_category_dict, fn_tcga_gtext_dict)

