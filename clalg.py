categories = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", 'root']


#The function evaluates level of taxonomic annotation for N% top hits for a gene to exclude bitscores 
#of low-level annotations and retain it in case of accidental absence of annotation

#For example "species", "genus", "family", "order", "class", "phylum", "superkingdom"
#               exist    exist       NA     exist    exist    exist     exist
#bitscore for NA will be accounted for
#For example "species", "genus", "family", "order", "class", "phylum", "superkingdom"
#               NA        NA       NA       exist     exist    exist     exist
#In this case bitscores for NA will be excluded

def evaluation(acc, index, bitscore, gene_dict, flag, taxonomy):
    if index < 8:
        if taxonomy[acc][categories[index]] == 'NA' and flag == 0:
            gene_dict[categories[index]] = 0
            index += 1
            return(evaluation(acc, index, bitscore, gene_dict, flag, taxonomy))
        elif taxonomy[acc][categories[index]] == 'NA' and flag == 1:
            gene_dict[categories[index]] = bitscore
            index += 1
            return(evaluation(acc, index, bitscore, gene_dict, flag, taxonomy))   
        else:
            gene_dict[categories[index]] = bitscore
            index += 1
            flag = 1 
            return(evaluation(acc, index, bitscore, gene_dict, flag, taxonomy))   
    else:
        return(gene_dict)


#Same function as previous, but for analysis of gene annotation after identification of LCA
def evaluation2(genes, index, gene_dict, flag):
    if index < 8:
        if genes[categories[index]] == 'NA' and flag == 0:
            gene_dict[categories[index]] = 0
            index += 1
            return(evaluation2(genes, index, gene_dict, flag))
        elif genes[categories[index]] == 'NA' and flag == 1:
            gene_dict[categories[index]] = 1
            index += 1
            return(evaluation2(genes, index, gene_dict, flag)) 
        else:
            gene_dict[categories[index]] = 1
            index += 1
            flag = 1 
            return(evaluation2(genes, index, gene_dict, flag))   
    else:
        return(gene_dict)

#Identification of LCA for N% of top hits
def calculation(eval_dict, contig, taxonomy):
    lca = ''   
    loc_dict = {}

    for cat in categories:
        if [taxonomy[item][cat] for item in eval_dict if eval_dict[item][cat] != 0]:
            l, num = zip(*[(taxonomy[item][cat], item) for item in eval_dict if eval_dict[item][cat] != 0])
            if len(set(l)) == 1 and l[0] != 'NA':
                lca  = cat
                b = [eval_dict[item][cat] for item in eval_dict if eval_dict[item][cat] != 0]
                bitscore = sum(b) / len(num)
                break
    
    for i in range(8):
        if i >= categories.index(lca):
            loc_dict[categories[i]] = taxonomy[num[0]][categories[i]]
        else:
            loc_dict[categories[i]] = 'NA'    
    loc_dict["bitscore"] = bitscore 
    return(loc_dict)

#Main function, contig classification for threshold >= 50%         
def classification_50(contig_dict, cont_bitscore):
    eval_dict = {}
    for gene in contig_dict:        
        eval_dict[gene] = evaluation2(contig_dict[gene], 0, {}, 0)
    loc_dict = {}
    loc_dict_2 = {}
    p = [contig_dict[item]['bitscore'] for item in contig_dict]
    B_threshold = sum(p) * cont_bitscore
    B_max = sum(p)
    for cat in categories:
        types = [contig_dict[item][cat] for item in contig_dict]    
        for name in set(types):
            l = [contig_dict[item]['bitscore'] for item in contig_dict if (contig_dict[item][cat] == name and eval_dict[item][cat] == 1)]
            if sum(l) >= B_threshold:             
                loc_dict[cat] = name + ':' + str(round(sum(l)/B_max*100, 1))
                loc_dict_2[cat] = name
    for cat in categories: 
        if not cat in loc_dict:
            loc_dict[cat] = "unclassified"
            loc_dict_2[cat] = "unclassified"      
    return(loc_dict, loc_dict_2) 

#Main function, contig classification for threshold < 50%                 
def classification_less(contig_dict, cont_bitscore):
    eval_dict = {}
    for gene in contig_dict:        
        eval_dict[gene] = evaluation2(contig_dict[gene], 0, {}, 0)
    loc_dict = {}
    p = [contig_dict[item]['bitscore'] for item in contig_dict]
    B_threshold = sum(p) * cont_bitscore
    B_max = sum(p)
    for cat in categories:
        types = [contig_dict[item][cat] for item in contig_dict]
        lu = []
        for name in set(types):
            l = [contig_dict[item]['bitscore'] for item in contig_dict if (contig_dict[item][cat] == name and eval_dict[item][cat] == 1)]
            if sum(l) >= B_threshold:
                p = name + ':' + str(round(sum(l)/B_max*100, 1))             
                lu.append(p)
        if len(lu) > 1:     
            loc_dict[cat] = ';'.join(lu)
        elif len(lu) == 1:
            loc_dict[cat] = str(lu).strip("['']")
    for cat in categories: 
        if not cat in loc_dict:
             loc_dict[cat] = "unclassified"      
    return(loc_dict) 
  
#Main function, contig classification for threshold >= 50%         
def nano_classification_50(contig_dict, B_threshold, num_fragments):
    eval_dict = {}
    for gene in contig_dict:        
        eval_dict[gene] = evaluation2(contig_dict[gene], 0, {}, 0)
    loc_dict = {}
    loc_dict_2 = {}
    for cat in categories[:-1]:
        types = [contig_dict[item][cat] for item in contig_dict]    
        for name in set(types):
            l = [contig_dict[item]['bitscore'] for item in contig_dict if (contig_dict[item][cat] == name and eval_dict[item][cat] == 1)]
            if len(l) / num_fragments > B_threshold:
                loc_dict[cat] = name + ':' + str(round(len(l)/num_fragments*100, 1))
                loc_dict_2[cat] = name
    loc_dict['root'] = '131567_cellular organisms:100.0'
    for cat in categories: 
        if not cat in loc_dict:
            loc_dict[cat] = "unclassified"
            loc_dict_2[cat] = "unclassified"      
    return(loc_dict, loc_dict_2) 

#Main function, contig classification for threshold < 50%                 
def nano_classification_less(contig_dict, B_threshold, num_fragments):
    eval_dict = {}
    for gene in contig_dict:        
        eval_dict[gene] = evaluation2(contig_dict[gene], 0, {}, 0)
    loc_dict = {}
    for cat in categories[:-1]:
        types = [contig_dict[item][cat] for item in contig_dict]
        lu = []
        for name in set(types):
            l = [contig_dict[item]['bitscore'] for item in contig_dict if (contig_dict[item][cat] == name and eval_dict[item][cat] == 1)]
            if len(l) / num_fragments >= B_threshold:
                p = name + ':' + str(round(len(l) / num_fragments*100, 1))             
                lu.append(p)
        if len(lu) > 1:     
            loc_dict[cat] = ';'.join(lu)
        elif len(lu) == 1:
            loc_dict[cat] = str(lu).strip("['']")
    loc_dict['root'] = '131567_cellular organisms:100.0'
    for cat in categories: 
        if not cat in loc_dict:
             loc_dict[cat] = "unclassified"      
    return(loc_dict) 
