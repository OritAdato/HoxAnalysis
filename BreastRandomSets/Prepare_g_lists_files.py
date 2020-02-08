import json
import os.path
import gzip

source_dir = "/home/HOX/Xena"
norm_cFromZena_fn = "TcgaTargetGtex_RSEM_Hugo_norm_count.gz"   #downloaded from Xena
TCGA_BRCA_exp_fromXena_fn = os.path.join(source_dir,"TCGA_BRCAexpression","HiSeqV2.gz")  # downloaded from Xena

gene_list_XenaPancanExp_fn = os.path.join(source_dir,"testRandom","gene_list_XenaPancanExp.json")
gene_list_XenaTCGAbrcaExp_fn = os.path.join(source_dir,"testRandom","gene_list_XenaTCGAbrcaExp.json")


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

norm_count_fn = os.path.join(source_dir,norm_cFromZena_fn)

gene_list_XenaPancanExp = prepare_gene_list(norm_count_fn)
write_l2_json(gene_list_XenaPancanExp, gene_list_XenaPancanExp_fn)

gene_list_XenaTCGAbrcaExp = prepare_gene_list(TCGA_BRCA_exp_fromXena_fn)
write_l2_json(gene_list_XenaTCGAbrcaExp, gene_list_XenaTCGAbrcaExp_fn)


