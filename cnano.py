import tax
import clalg

categories2 = ["superkingdom","phylum","class","order","family","genus","species"]
categories = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", 'root']

def cat_without_orfs(path_to_taxonomy_files, fasta, prefix, fragment, diamond_file, cust_bit_thresh, qual_threshold, cont_bitscore, no_summary, contigs_info, length_dict, log, quiet):

    main_dict, taxonomy = tax.parsing_alignment_long(path_to_taxonomy_files, diamond_file, cust_bit_thresh, qual_threshold, log, quiet)

    file_out = open(prefix + '_contigs_classification.txt', 'w')
    if cont_bitscore >= 0.5:
        file_out.write('#Contig_ID\tClassification_level\tClassification\tNo_annotated_fragments\tNo_fragments\tAnnontated_proportion\ttaxid_superkingdom:percentage\ttaxid_phylum:percentage\ttaxid_class:percentage\ttaxid_order:percentage\ttaxid_family:percentage\ttaxid_genus:percentage\ttaxid_species:percentage' + '\n')
    else:
        file_out.write('#Contig_ID\tNo_annotated_fragments\tNo_fragments\tAnnontated_proportion\ttaxid_superkingdom:percentage\ttaxid_phylum:percentage\ttaxid_class:percentage\ttaxid_order:percentage\ttaxid_family:percentage\ttaxid_genus:percentage\ttaxid_species:percentage' + '\n')

    contigs_names2 = sorted([i for i in contigs_info])
    summary = {}
    for cat in categories:
        summary[cat] = {}
        
    if quiet == False:
        print('CAT starts classification of contigs...')
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
            #Contig classification
            if cont_bitscore >= 0.5 and contig_dict:
                contig_annot, contig_annot2 = clalg.nano_classification_50(contig_dict, cont_bitscore, length_dict[contig])
                output = []
                output.append(contig)
                for cat in categories:
                    if contig_annot[cat] != 'unclassified':
                        output.append(cat)
                        output.append(contig_annot[cat])
                        break
                output.append(str(number_genes_annotated))
                output.append(str(length_dict[contig]))
                output.append(str(round(number_genes_annotated/length_dict[contig] * 100, 2)))
                for cat in categories2:
                    if contig_annot2[cat] in summary[cat]:
                        summary[cat][contig_annot2[cat]][0] += 1
                        summary[cat][contig_annot2[cat]][1] += contigs_info[contig]['length']
                    else:
                        summary[cat][contig_annot2[cat]] = [0, 0]
                        summary[cat][contig_annot2[cat]][0] = 1
                        summary[cat][contig_annot2[cat]][1] += contigs_info[contig]['length']
                    output.append(contig_annot[cat])
                file_out.write('\t'.join(output) + '\n')
            elif cont_bitscore < 0.5 and contig_dict:
                contig_annot = clalg.nano_classification_less(contig_dict, cont_bitscore, length_dict[contig])
                output = []
                output.append(contig)
    
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
            file7.write('#Total number of annotated contigs is ' + str(len(main_dict)) + '\n')
            file7.write('#Total number of unannotated contigs is ' + str(len(contigs_info) - len(main_dict)) + '\n')
            file7.write('#Classification_level' + '\t' + 'Taxon_name' + '\t' + 'Number_of_contigs' + '\t' + 'Sum_length_of_contigs' + '\n')
            for cat in categories2:  
                sum_tot = []      
                for item in summary[cat]:
                    l = [cat, item, summary[cat][item][0], summary[cat][item][1]]
                    sum_tot.append(l)
                total = sorted(sum_tot, key=lambda x: x[2], reverse=True)
                for item in total:
                    file7.write('\t'.join([str(x) for x in item]) + '\n')

