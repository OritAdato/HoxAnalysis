def prep_sample_category_dict(category_tab_fn):
    samp_category_dict = {}
    with open(category_tab_fn,'r') as category_tab_f:
        for line in category_tab_f:
            line_l = line.split('\t')
            samp = line_l[0]
            samp_category = line_l[1].rstrip('\n\r')
            samp_category_dict.setdefault(samp,[]).append(samp_category)

    return samp_category_dict

def prep_list_of_tumorHealth_TCGA_GTEX(tcga_allBreast_rownames):

    TCGA_breast_tumor_sampl = []
    GTEX_breast_health_sampl = []
    for samp in tcga_allBreast_rownames:
        sampl = samp.strip('\r\n')
        if sampl[:4] == 'TCGA':
            TCGA_breast_tumor_sampl.append(sampl)
        else:
            if sampl[:4] == 'GTEX':
                GTEX_breast_health_sampl.append(sampl)
            else:
                print " ====> Error Error - unexpected sample type: ", sampl
                exit()

    return TCGA_breast_tumor_sampl, GTEX_breast_health_sampl

def prep_list_of_tumorHealth_TCGAonly(tcga_allBreast_rownames):

    TCGA_breast_tumor_sampl = []
    TCGA_breast_health_sampl = []
    for samp in tcga_allBreast_rownames:
        sampl = samp.strip('\r\n')
        if sampl[len(sampl)-3:len(sampl)] != '-11':
            TCGA_breast_tumor_sampl.append(sampl)
        else:
            TCGA_breast_health_sampl.append(sampl)

    TCGA_breast_tumor_samplNoH=remove_health_samplesFrom_list(TCGA_breast_tumor_sampl)
    return TCGA_breast_tumor_samplNoH, TCGA_breast_health_sampl

def prepare_gene_list(norm_count_fn):
#
# gene name is first column in file
#
    g_list =[]
    with gzip.open(norm_count_fn, 'rb') as norm_count_f:
        ln = 0
        for gline in norm_count_f:
            line_l = gline.split('\t')
            if line_l[0] == 'sample':
                continue
            else:
                gname = line_l[0]
                g_list.append(gname)
    return g_list

def write_l2_json(g_list, gene_list_XenaPancanExp_fn):
    with open(gene_list_XenaPancanExp_fn,'w') as gene_list_XenaPancanExp_f:
        json.dump(g_list,gene_list_XenaPancanExp_f)
    return

def join2lists(TCGA_breast_tumor_100_randSampl,GTEX_breast_health_100_randSampl):
    oneBiglist = []
    map(oneBiglist.extend, [TCGA_breast_tumor_100_randSampl,GTEX_breast_health_100_randSampl])
    return oneBiglist

def remove_health_samplesFrom_list(TCGA_sample_list):
    tumor_list = []
    for samp_name in TCGA_sample_list:
        if samp_name[len(samp_name)-3:len(samp_name)] != '-11':
            tumor_list.append(samp_name)
    return tumor_list

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
