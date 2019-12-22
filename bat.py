import subprocess
import clalg
import tax
import json
import copy

categories2 = ["superkingdom","phylum","class","order","family","genus","species"]
categories = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", 'root']


def bat_with_orfs(path_to_taxonomy_files, prefix, diamond_file, protein_fasta, qual_threshold, cust_bit_thresh, cont_bitscore, no_lca, no_summary, prot_con_links, contigs_info, prot_weight, bins_dict, log, quiet, bins_b2_conf, bins_b2_other):
    #Calculate contribution of annotated ORFs to taxonomic classification 
    def weights_comp(contig, contig_dict):
        Art_Bmax = 0
        Art_Bann = 0
        for scaf in prot_weight[contig]:
            Art_Bmax += prot_weight[contig][scaf]
        for genes in contig_dict:
            Art_Bann += prot_weight[contig][genes]
        return(str(round(Art_Bann / Art_Bmax * 100, 2)))

    def weights_comp_bins(contig, contig_dict):
        Art_Bmax = 0
        Art_Bann = 0
        for cont in bins_dict[contig]:
            for scaf in prot_weight[cont]:
                Art_Bmax += prot_weight[cont][scaf]
                if scaf in contig_dict:
                    Art_Bann += prot_weight[cont][scaf]
        return(str(round(Art_Bann / Art_Bmax * 100, 2)))

    main_dict, taxonomy = tax.parsing_alignment(path_to_taxonomy_files, diamond_file, cust_bit_thresh, qual_threshold, prot_con_links, log, quiet)
    comm = ('mkdir ' + prefix + '_contigs_classification_per_bin').split()
    subprocess.call(comm)
    
    if no_lca == False:
        subprocess.call(['echo', '-e', '#Contig_ID\tAverage_bitscore\ttaxid_superkingdom\ttaxid_phylum\ttaxid_class\ttaxid_order\ttaxid_family\ttaxid_genus\ttaxid_species'], stdout=open(prefix + '_LCA_per_ORF.txt','w'))


    summary = {}
    for cat in categories:
        summary[cat] = {}


    if quiet == False:
        print('CAT starts classification of contigs individually')
    log.write('CAT starts classification of contigs individually' + '\n')

    for file_name in sorted(bins_dict):
        file_out = open('./' + prefix + '_contigs_classification_per_bin/' + file_name + '_contigs_classification.txt', 'w')
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

        for contig in sorted(bins_dict[file_name]):
            if contig in main_dict:
                gene_dict = copy.deepcopy(main_dict[contig])
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

                    output.append(str(len(bins_dict[file_name][contig])))
                    if protein_fasta:
                        weight = weights_comp(contig, contig_dict)
                        output.append(weight)
                    output.append(str(number_genes_annotated))
                    for cat in categories2:
                        if contig_annot2[cat] in summary[cat]:
                            summary[cat][contig_annot2[cat]][0] += 1
                            summary[cat][contig_annot2[cat]][1] += len(bins_dict[file_name][contig])
                            summary[cat][contig_annot2[cat]][2] += contigs_info[contig]['length']
                        else:
                            summary[cat][contig_annot2[cat]] = [0, 0, 0]
                            summary[cat][contig_annot2[cat]][0] = 1
                            summary[cat][contig_annot2[cat]][1] += len(bins_dict[file_name][contig])
                            summary[cat][contig_annot2[cat]][2] += contigs_info[contig]['length']
                        output.append(contig_annot[cat])
                    file_out.write('\t'.join(output) + '\n')
                elif cont_bitscore < 0.5 and contig_dict:
                    contig_annot = clalg.classification_less(contig_dict, cont_bitscore)
                    output = []
                    output.append(contig)
                    output.append(str(len(bins_dict[file_name][contig])))
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

        sum_file = open(prefix + '_contigs_summary.txt', 'w')
    #Generation of summary file
    if cont_bitscore >= 0.5 and not no_summary:
        sum_file.write('#Total number of contigs is ' + str(len(contigs_info)) + '\n')
        sum_file.write('#Number of annotated contigs is ' + str(len(main_dict)) + '\n')
        sum_file.write('#Number of unannotated contigs is ' + str(len(contigs_info) - len(main_dict)) + '\n')
        sum_file.write('#Classification_level' + '\t' + 'Taxon_name' + '\t' + 'Number_of_contigs' + '\t' + 'Sum_number_of_proteins' + '\t' + 'Sum_length_of_contigs' + '\n')
        for cat in categories2:  
            sum_tot = []      
            for item in summary[cat]:
                l = [cat, item, summary[cat][item][0], summary[cat][item][1], summary[cat][item][2]]
                sum_tot.append(l)
            total = sorted(sum_tot, key=lambda x: x[2], reverse=True)
            for item in total:
                sum_file.write('\t'.join([str(x) for x in item]) + '\n')
    sum_file.close()
    
    
    if quiet == False:
        print('CAT is starting classification of bins...')
    log.write('CAT is starting classification of bins...\n')
    
     #### Classification of bin as a whole
    no_hits = 0
    summary = {}
    for cat in categories:
        summary[cat] = {}

    file_out = open(prefix + '_bins_classification.txt', 'w')
    if protein_fasta:
        file_out.write('#Bin_ID\tClassification_level\tClassification\tNo_of_contigs_in_bin\tNo_genes\tContribution_of_annotated_ORFs\tNo_annotated_genes\tConfident_taxid_superkingdom:contribution\tConfident_taxid_phylum:contribution\tConfident_taxid_class:contribution\tConfident_taxid_order:contribution\tConfident_taxid_family:contribution\tConfident_taxid_genus:contribution\tConfident_taxid_species:contribution\tOther_taxid_superkingdom:contribution\tOther_taxid_phylum:contribution\tOther_taxid_class:contribution\tOther_taxid_order:contribution\tOther_taxid_family:contribution\tOther_taxid_genus:contribution\tOther_taxid_species:contribution' + '\n')
    else:
        file_out.write('#Bin_ID\tClassification_level\tClassification\tNo_of_contigs_in_bin\tNo_genes\tNo_annotated_genes\tConfident_taxid_superkingdom:contribution\tConfident_taxid_phylum:contribution\tConfident_taxid_class:contribution\tConfident_taxid_order:contribution\tConfident_taxid_family:contribution\tConfident_taxid_genus:contribution\tConfident_taxid_species:contribution\tOther_taxid_superkingdom:contribution\tOther_taxid_phylum:contribution\tOther_taxid_class:contribution\tOther_taxid_order:contribution\tOther_taxid_family:contribution\tOther_taxid_genus:contribution\tOther_taxid_species:contribution' + '\n')

    for contig in sorted(bins_dict):
        cont_bitscore = bins_b2_conf
        gene_dict = {}
        for contig2 in bins_dict[contig]:
            if contig2 in main_dict:
                for gene in main_dict[contig2]:
                    gene_dict[gene] = copy.deepcopy(main_dict[contig2][gene])
        if gene_dict:
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
        #Contig classification
            if contig_dict:
                contig_annot, contig_annot2 = clalg.classification_50(contig_dict, cont_bitscore)
                output = []
                output.append(contig)
                for cat in categories:
                    if contig_annot[cat] != 'unclassified':
                        output.append(cat)
                        output.append(contig_annot[cat])
                        break
                output.append(str(len(bins_dict[contig])))
                output.append(str(len([bins_dict[contig][contig2][gene] for contig2 in bins_dict[contig] for gene in bins_dict[contig][contig2]])))
                if protein_fasta:
                    weight = weights_comp_bins(contig, contig_dict)
                    output.append(weight)
                output.append(str(number_genes_annotated))
                for cat in categories2:
                    if contig_annot2[cat] in summary[cat]:
                        summary[cat][contig_annot2[cat]][0] += 1
                        for true_contig in bins_dict[contig]:
                            summary[cat][contig_annot2[cat]][2] += contigs_info[true_contig]['length']
                            summary[cat][contig_annot2[cat]][1] += len(bins_dict[contig][true_contig])
                    else:
                        summary[cat][contig_annot2[cat]] = [0, 0, 0]
                        summary[cat][contig_annot2[cat]][0] = 1
                        for true_contig in bins_dict[contig]:
                            summary[cat][contig_annot2[cat]][1] += len(bins_dict[contig][true_contig])
                            summary[cat][contig_annot2[cat]][2] += contigs_info[true_contig]['length']
                    output.append(contig_annot[cat])
                cont_bitscore = bins_b2_other
                contig_annot = clalg.classification_less(contig_dict, cont_bitscore)
                for cat in categories2:
                    output.append(contig_annot[cat])
                file_out.write('\t'.join(output) + '\n')
        else:
            file_out.write(contig + '\tNo_hits')
            no_hits += 1
    file_out.close()
    
    sum_file = open(prefix + '_bins_summary.txt', 'w')
#Generation of summary file

    if not no_summary:
        sum_file.write('#Total number of bins is ' + str(len(bins_dict)) + '\n')
        sum_file.write('#Number of annotated bins is ' + str(len(bins_dict) - no_hits) + '\n')
        sum_file.write('#Number of unannotated bins is ' + str(no_hits) + '\n')
        sum_file.write('#Classification_level' + '\t' + 'Taxon_name' + '\t' + 'Number_of_bins' + '\t' + 'Sum_number_of_proteins' + '\t' + 'Sum_length_of_contigs_in_bins' + '\n')
        for cat in categories2:  
            sum_tot = []      
            for item in summary[cat]:
                l = [cat, item, summary[cat][item][0], summary[cat][item][1], summary[cat][item][2]]
                sum_tot.append(l)
            total = sorted(sum_tot, key=lambda x: x[2], reverse=True)
            for item in total:
                sum_file.write('\t'.join([str(x) for x in item]) + '\n')
    sum_file.close()
