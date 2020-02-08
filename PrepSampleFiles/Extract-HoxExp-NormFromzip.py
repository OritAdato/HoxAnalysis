'''
Created on Nov 5, 2017

@author: orit
'''
import gzip
import os.path
import pprint
#  Extract only HOX genes

def prepare_hox_glist(hox_glist_fn):
    title = True
    hox_list =[]
    AS_hox = []
    with open(hox_glist_fn,'r') as hox_glist_f:
        for line in hox_glist_f:
            if title:
                title = False            
                continue
            hox_line = line.split(',')
#            hox_list.append(hox_line[0])
            if '-' in hox_line[0]:    # this condition filters the AS Hoxes
                AS_hox.append(hox_line[0])
            else:
                hox_list.append(hox_line[0])  
    for AS_hox_gene in AS_hox:
        hox_list.append(AS_hox_gene)             # AS will be appended at the end of the list 
    
    return hox_list

def prep_sample_category_dict(category_tab_fn):
    samp_category_dict = {}
    with open(category_tab_fn,'r') as category_tab_f:
        for line in category_tab_f:
            line_l = line.split('\t')
            samp = line_l[0]
            samp_category = line_l[1].rstrip('\n\r')
            samp_category_dict.setdefault(samp,[]).append(samp_category)            

    return samp_category_dict

def print_title_line(tcga_gtex_exp_f,hox_genes_l):
    title_line = 'Samp ID, Cohort,Category'
    for ind in range(len(hox_genes_l)):
        title_line = title_line + ',' + hox_genes_l[ind]
    title_line = title_line + '\n'
    tcga_gtex_exp_f.write(title_line)
    return

def prepare_line4print (samp,samp_category,gene_exp_l):
    dashpos = samp.find('-')
    cohort = samp[0:dashpos]
    line4p = samp + ',' + str(cohort) + ',' + str(samp_category)
    for ind in range(len(gene_exp_l)):
#        print 'ind:', str(ind)
        gene_exp = gene_exp_l[ind][0]
        line4p = line4p + ',' + str(gene_exp)
    line4p = line4p + '\n'
    return line4p

def read_file_det(norm_count_fn,tcga_gtex_exp_f,hox_genes_l):
    
    with gzip.open(norm_count_fn, 'rb') as norm_count_f:
        samp_list =[]
        hox_exp_dict = {}
        for gline in norm_count_f:
            line_l = gline.split('\t')
            if line_l[0] == 'sample':
                for ind in range (1, len(line_l)):
                    samp_id = line_l[ind].strip()
                    samp_list.append(samp_id)
            else:
                gname = line_l[0]
                if gname in hox_genes_l:
                    for ind in range(1,len(line_l)):
                        g_exp= line_l[ind].strip()
                        hox_exp_dict.setdefault(gname,[]).append(gexp)
 #       for hox_g in hox_genes_l:
                    
    return

def read_file(norm_count_fn,tcga_gtex_exp_f,hox_genes_l):
    
    with gzip.open(norm_count_fn, 'rb') as norm_count_f:
        for gline in norm_count_f:
            line_l = gline.split('\t')
            if line_l[0] == 'sample':
                tcga_gtex_exp_f.write(gline)
            else:
                gname = line_l[0]
                if gname in hox_genes_l:
                    tcga_gtex_exp_f.write(gline)      
                    
    return

hox_glist_fn = "/old_data/orit/HOX/Hox-Identifiers4Xena.csv"   # For HOX genes

category_tab_fn = "/data/old_data/orit/HOX/Xena/TCGA_GTEX_category.txt/TCGA_GTEX_category.txt"
tcga_gtex_exp_fn = "/data/old_data/orit/HOX/Xena/TCGA_GTEX_exp_norm_count.csv"  # HOX genes
source_dir = "/data/old_data/orit/HOX/Xena"
norm_cFromZena_fn = "TcgaTargetGtex_RSEM_Hugo_norm_count.gz"

hox_genes_l = prepare_hox_glist(hox_glist_fn)
samp_category_dict =  prep_sample_category_dict(category_tab_fn)
norm_count_fn = os.path.join(source_dir,norm_cFromZena_fn)

tcga_gtex_exp_f = open(tcga_gtex_exp_fn,'w')
#title_line = print_title_line(tcga_gtex_exp_f,hox_genes_l)
read_file(norm_count_fn,tcga_gtex_exp_f,hox_genes_l) 
tcga_gtex_exp_f.close()
exit()


