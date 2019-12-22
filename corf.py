import subprocess
import clalg
import tax
import json

categories2 = ["superkingdom","phylum","class","order","family","genus","species"]
categories = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", 'root']


def cat_with_orfs(path_to_taxonomy_files, fasta, prefix, diamond_file, protein_fasta, qual_threshold, cust_bit_thresh, cont_bitscore, no_lca, no_summary, prot_con_links, contigs_info, prot_weight, con_prot_links, log, quiet):
    #Calculate contribution of annotated ORFs to taxonomic classification 
    def weights_comp(contig, contig_dict):
        Art_Bmax = 0
        Art_Bann = 0
        for scaf in prot_weight[contig]:
            Art_Bmax += prot_weight[contig][scaf]
        for genes in contig_dict:
            Art_Bann += prot_weight[contig][genes]
        return(str(round(Art_Bann / Art_Bmax * 100, 2)))
    main_dict, taxonomy = tax.parsing_alignment(path_to_taxonomy_files, diamond_file, cust_bit_thresh, qual_threshold, prot_con_links, log, quiet)

    if no_lca == False:
        subprocess.call(['echo', '-e', '#Contig_ID\tAverage_bitscore\ttaxid_superkingdom\ttaxid_phylum\ttaxid_class\ttaxid_order\ttaxid_family\ttaxid_genus\ttaxid_species'], stdout=open(prefix + '_LCA_per_ORF.txt','w'))
                         
    #Calculate contribution of annotated ORFs to taxonomic classification 

    contigs_names2 = sorted([i for i in contigs_info])
    file_out = open(prefix + '_contigs_classification.txt', 'w') 
    
    if protein_fasta:
        if cont_bitscore >= 0.5:
            file_out.write('#Contig_ID\tClassification_level\tClassification\tNo_genes\tContribution_of_annotated_ORFs\tNo_annotated_genes\ttaxid_superkingdom:contribution\ttaxid_phylum:contribution\ttaxid_class:contribution\ttaxid_order:contribution\ttaxid_family:contribution\ttaxid_genus:contribution\ttaxid_species:contribution' + '\n')
        else:
            file_out.write('#Contig_ID\tNo_genes\tContribution_of_annotated_ORFs\tNo_annotated_genes\ttaxid_superkingdom:contribution\ttaxid_phylum:contribution\ttaxid_class:contribution\ttaxid_order:contribution\ttaxid_family:contribution\ttaxid_genus:contribution\ttaxid_species:contribution' + '\n')
    else:
        if cont_bitscore >= 0.5:
            file_out.write('#Contig_ID\tClassification_level\tClassification\tNo_genes\tNo_annotated_genes\ttaxid_superkingdom:contribution\ttaxid_phylum:contribution\ttaxid_class:contribution\ttaxid_order:contribution\ttaxid_family:contribution\ttaxid_genus:contribution\ttaxid_species:contribution' + '\n')
        else:
            file_out.write('#Contig_ID\tNo_genes\tNo_annotated_genes\ttaxid_superkingdom:contribution\ttaxid_phylum:contribution\ttaxid_class:contribution\ttaxid_order:contribution\ttaxid_family:contribution\ttaxid_genus:contribution\ttaxid_species:contribution' + '\n')
    
    summary = {}
    for cat in categories:
        summary[cat] = {}
                   		            
    if quiet == False:
        print('CAT starts classification of contigs')
    log.write('CAT starts classification of contigs' + '\n')

    for contig in contigs_names2:
        if contig in main_dict:
            gene_dict = main_dict[contig]
            contig_dict = {} 
        #Evaluation of annotations per gene and analysing taxonomic 'level' of annotation       
            gene_dict_eval = {}
            for gene in gene_dict:
                gene_dict_eval[gene] = {}            
                if len(gene_dict[gene]['accessions']) > 1:
                    for item in range(len(gene_dict[gene]['accessions'])):
                        if taxonomy[gene_dict[gene]['accessions'][item]]:
                            gene_dict_eval[gene][gene_dict[gene]['accessions'][item]] = clalg.evaluation(gene_dict[gene]['accessions'][item], int(0), gene_dict[gene]['bitscores'][item], {}, 0, taxonomy)
                    if not gene_dict_eval[gene]:
                        gene_dict_eval.pop(gene, None)
                    else:        
                        contig_dict[gene] = clalg.calculation(gene_dict_eval[gene], contig, taxonomy)
                elif len(gene_dict[gene]['accessions']) == 1 and taxonomy[gene_dict[gene]['accessions'][0]]:                                   
                    gene_dict_eval[gene][gene_dict[gene]['accessions'][0]] = clalg.evaluation(gene_dict[gene]['accessions'][0], int(0), gene_dict[gene]['bitscores'][0], {}, 0, taxonomy)       
                    loc_dict = {}
                    for cat in categories:           
                        if gene_dict_eval[gene][gene_dict[gene]['accessions'][0]][cat] != 0:
                            loc_dict[cat] = taxonomy[gene_dict[gene]['accessions'][0]][cat]                       
                        else:
                            loc_dict[cat] = "NA"
                    loc_dict['bitscore'] = gene_dict_eval[gene][gene_dict[gene]['accessions'][0]]['root']
                    contig_dict[gene] = loc_dict       
            number_genes_annotated = len(contig_dict)
            if number_genes_annotated > 0 and not no_lca:
                out_orf = open(prefix + '_LCA_per_ORF.txt', 'a')
                for contig_orf in sorted(contig_dict):
                    loc_list = []
                    loc_flag = 0
                    for cat in categories[:-1]:
                        if contig_dict[contig_orf][cat] == 'NA' and loc_flag == 0:
                            loc_list.append('unclassified')
                        elif contig_dict[contig_orf][cat] == 'NA' and loc_flag == 1:
                            loc_list.append('NA')
                        else:
                            loc_list.append(contig_dict[contig_orf][cat])
                            loc_flag = 1
                    loc_list = loc_list[::-1]
                    loc_list.insert(0, str(round(contig_dict[contig_orf]['bitscore'], 2)))
                    loc_list.insert(0, contig_orf)
                    out_orf.write('\t'.join(loc_list) + '\n')
                out_orf.close()
        
        #Contig classification
            if cont_bitscore >= 0.5 and contig_dict:
                contig_annot, contig_annot2 = clalg.classification_50(contig_dict, cont_bitscore)
                output = []
                output.append(contig)
                for cat in categories:
                    if contig_annot[cat] != 'unclassified':
                        output.append(cat)
                        output.append(contig_annot[cat])
                        break

                output.append(str(con_prot_links[contig]))
                if protein_fasta:
                    weight = weights_comp(contig, contig_dict)
                    output.append(weight)
                output.append(str(number_genes_annotated))
                for cat in categories2:
                    if contig_annot2[cat] in summary[cat]:
                        summary[cat][contig_annot2[cat]][0] += 1
                        summary[cat][contig_annot2[cat]][1] += con_prot_links[contig]
                        summary[cat][contig_annot2[cat]][2] += contigs_info[contig]['length']
                    else:
                        summary[cat][contig_annot2[cat]] = [0, 0, 0]
                        summary[cat][contig_annot2[cat]][0] = 1
                        summary[cat][contig_annot2[cat]][1] += con_prot_links[contig]
                        summary[cat][contig_annot2[cat]][2] += contigs_info[contig]['length']
                    output.append(contig_annot[cat])
                file_out.write('\t'.join(output) + '\n')
            elif cont_bitscore < 0.5 and contig_dict:
                contig_annot = clalg.classification_less(contig_dict, cont_bitscore)
                output = []
                output.append(contig)
                output.append(str(con_prot_links[contig]))
                if protein_fasta:
                    weight = weights_comp(contig, contig_dict)
                    output.append(weight)
                output.append(str(number_genes_annotated))
                for cat in categories2:
                    output.append(contig_annot[cat])
                file_out.write('\t'.join(output) + '\n')
            else:
                file_out.write(contig + '\tNo_hits' + '\n')
        else:
            file_out.write(contig + '\tNo_hits' + '\n') 
    file_out.close()

    #Generation of summary file

    if cont_bitscore >= 0.5 and not no_summary:
        with open(prefix + '_summary.txt', 'w') as file7:
            file7.write('#Total number of contigs is ' + str(len(contigs_info)) + '\n')
            file7.write('#Number of annotated contigs is ' + str(len(main_dict)) + '\n')
            file7.write('#Number of unannotated contigs is ' + str(len(contigs_info) - len(main_dict)) + '\n')
            file7.write('#Classification_level' + '\t' + 'Taxon_name' + '\t' + 'Number_of_contigs' + '\t' + 'Sum_number_of_proteins' + '\t' + 'Sum_length_of_contigs' + '\n')
            for cat in categories2:  
                sum_tot = []      
                for item in summary[cat]:
                    l = [cat, item, summary[cat][item][0], summary[cat][item][1], summary[cat][item][2]]
                    sum_tot.append(l)
                total = sorted(sum_tot, key=lambda x: x[2], reverse=True)
                for item in total:
                    file7.write('\t'.join([str(x) for x in item]) + '\n')

