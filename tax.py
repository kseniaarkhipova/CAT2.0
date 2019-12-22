import csv
import re
import json

categories = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", 'root']

def parsing_taxonomy_files(taxonomy, path_to_taxonomy_files, log, quiet):
    names_dict = {}
    rank_dict = {}
    parent_dict = {}
    def parse_nodes(taxid, data_dict):
        if taxid in parent_dict and taxid in rank_dict:
            while parent_dict[taxid] != taxid:
                data_dict[rank_dict[taxid]] = taxid + '_' + names_dict[taxid]
                taxid = parent_dict[taxid]
        else:
              data_dict = {}  
        return(data_dict)

    with open(path_to_taxonomy_files + 'names.dmp') as names_file:
        names = csv.reader(names_file,delimiter='\t',quoting=csv.QUOTE_NONE)     
        for line in names:
            if line[6] == 'scientific name':
                names_dict[line[0]] = line[2]

    with open(path_to_taxonomy_files + 'nodes.dmp','r')as nodes_file:
        nodes = csv.reader(nodes_file,delimiter='\t',quoting=csv.QUOTE_NONE)
        for line in nodes:
            parent_dict[line[0]] = line[2]
            rank_dict[line[0]] = line[4]
            
    if quiet == False:
        print('CAT is parsing taxonomy file with connections of accession identifiers to taxonomical identifiers. This will take approximally 5 minutes...')
    log.write('CAT is parsing taxonomy file with connections of accession identifiers to taxonomical identifiers. This will take approximally 5 minutes...\n')


    with open(path_to_taxonomy_files + 'prot.accession2taxid') as file1:
        acc_tax_csv=csv.reader(file1,delimiter='\t',quoting=csv.QUOTE_NONE)
        for line in acc_tax_csv:        
            if line[0] != 'accession':
                if line[0] in taxonomy:
                    z = parse_nodes(line[2], {})
                    if z:
                        for cat in categories:
                            if cat in z.keys():
                                taxonomy[line[0]][cat] = z[cat] 
                            else:
                                taxonomy[line[0]][cat] = 'NA'
                    else:
                        for cat in categories:
                            taxonomy[line[0]][cat] = 'NA'
                    taxonomy[line[0]]['root'] = '131567_cellular organisms'
    return(taxonomy)
                    
                    
                    
def parsing_alignment(path_to_taxonomy_files, diamond_file, cust_bit_thresh, qual_threshold, prot_con_links, log, quiet):
    #Parsing of alignment file
    if quiet == False:
        print('CAT is parsing file with alignment. This can take some time.')
    log.write('CAT is parsing file with alignment. This can take some time.\n')

    main_dict = {}
    with open(diamond_file) as file_1:
        alignment=csv.reader(file_1,delimiter='\t')
        taxonomy = {}
        for line in alignment:
            without_version = ''
            without_version = re.search('^(.+?)(?=(\.|:|$))', line[1]).group(1)
            if line[0] in prot_con_links:
                contig_name = prot_con_links[line[0]]
                gene_name = line[0]
                if without_version not in taxonomy:
                    taxonomy[without_version] = {}
                if contig_name not in main_dict:
                    main_dict[contig_name] = {}   
                if gene_name not in main_dict[contig_name] and float(line[-1]) > qual_threshold:
                    main_dict[contig_name][gene_name] = {}
                    main_dict[contig_name][gene_name]['accessions'] = [without_version]
                    main_dict[contig_name][gene_name]['bitscores'] = [float(line[-1])]
                    bit_threshold = float(line[-1]) * cust_bit_thresh
                elif gene_name in main_dict[contig_name] and float(line[-1]) > bit_threshold and float(line[-1]) > qual_threshold:            
                    main_dict[contig_name][gene_name]['accessions'].append(without_version)
                    main_dict[contig_name][gene_name]['bitscores'].append(float(line[-1]))
    taxonomy = parsing_taxonomy_files(taxonomy, path_to_taxonomy_files, log, quiet)
    return(main_dict, taxonomy)


def parsing_alignment_long(path_to_taxonomy_files, diamond_file, cust_bit_thresh, qual_threshold, log, quiet):
    def gene_chooser(dict_per_fragment):
        if len(dict_per_fragment) == 1:
            for i in dict_per_fragment:
                protein = i
        else:
            a = []
            b = []
            for frame in dict_per_fragment:
                a.append(frame)
                b.append(dict_per_fragment[frame]['bitscore'][0])
            protein = a[b.index(max(b))]
        bit_threshold = dict_per_fragment[protein]['bitscore'][0] * cust_bit_thresh
        new_dict = {}
        new_dict['acc'] = []
        new_dict['bitscore'] = []
        for i in range(len(dict_per_fragment[protein]['bitscore'])):
            if dict_per_fragment[protein]['bitscore'][i] >= bit_threshold:
                new_dict['acc'].append(dict_per_fragment[protein]['acc'][i])
                new_dict['bitscore'].append(dict_per_fragment[protein]['bitscore'][i])
        return(protein, new_dict)

    #Parsing of alignment file
    if quiet == False:
        print('CAT is parsing file with alignment. This can take some time.')
    log.write('CAT is parsing file with alignment. This can take some time.\n')

    
    main_dict = {}
    with open(diamond_file) as file_1:
        alignment=csv.reader(file_1,delimiter='\t')
        taxonomy = {}
        dict_per_fragment = {}
        previous_gene = ''
        previous_contig = ''
        for line in alignment:
            contig_name, fragment, frame = re.search('(\S+)_(fragment_\d+)_(frame_\d)', line[0]).group(1,2,3)
            gene_name = contig_name + '_' + fragment
            without_version = re.search('^(.+?)(?=(\.|:|$))', line[1]).group(1)
            if not previous_gene and float(line[-1]) > qual_threshold:
                previous_gene = gene_name
                previous_contig = contig_name
                if frame not in dict_per_fragment:
                    dict_per_fragment[frame] = {}
                    dict_per_fragment[frame]['acc'] = []
                    dict_per_fragment[frame]['bitscore'] = []
                dict_per_fragment[frame]['acc'].append(without_version)
                dict_per_fragment[frame]['bitscore'].append(float(line[-1]))
            elif previous_gene and previous_gene == gene_name and float(line[-1]) > qual_threshold:
                if frame not in dict_per_fragment:
                    dict_per_fragment[frame] = {}
                    dict_per_fragment[frame]['acc'] = []
                    dict_per_fragment[frame]['bitscore'] = []

                dict_per_fragment[frame]['acc'].append(without_version)
                dict_per_fragment[frame]['bitscore'].append(float(line[-1]))
            elif previous_gene and previous_gene != gene_name and float(line[-1]) > qual_threshold:
                frame_chosen,new_dict = gene_chooser(dict_per_fragment)
                for acc in new_dict['acc']:
                    taxonomy[acc] = {}
                gene_chosen = previous_gene + '_' + frame_chosen
                if previous_contig not in main_dict:
                    main_dict[previous_contig] = {}
                main_dict[previous_contig][gene_chosen] = {}
                main_dict[previous_contig][gene_chosen]['accessions'] = new_dict['acc']
                main_dict[previous_contig][gene_chosen]['bitscores'] = new_dict['bitscore']
                dict_per_fragment = {}
                dict_per_fragment[frame] = {}
                dict_per_fragment[frame]['acc'] = []
                dict_per_fragment[frame]['bitscore'] = []
                dict_per_fragment[frame]['acc'].append(without_version)
                dict_per_fragment[frame]['bitscore'].append(float(line[-1]))
                previous_gene = gene_name
                previous_contig = contig_name

    frame_chosen, new_dict = gene_chooser(dict_per_fragment)
    for acc in new_dict['acc']:
        taxonomy[acc] = {}
    gene_chosen = previous_gene + '_' + frame_chosen
    if previous_contig not in main_dict:
        main_dict[previous_contig] = {}
    main_dict[previous_contig][gene_chosen] = {}
    main_dict[previous_contig][gene_chosen]['accessions'] = new_dict['acc']
    main_dict[previous_contig][gene_chosen]['bitscores'] = new_dict['bitscore']
    taxonomy = parsing_taxonomy_files(taxonomy, path_to_taxonomy_files, log, quiet)
    return(main_dict, taxonomy)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
