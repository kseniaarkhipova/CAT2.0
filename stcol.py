import os
import re

sub = {'A' : 4, 'R' : 5, 'D' : 6, 'N' : 6, 'C' : 9, 'E' : 5, 'Q' : 5, 'G' : 6, 'H' : 8, 'I' : 4, 'L' : 4, 'K' : 5, 'M' : 5, 'F' : 6, 'P' : 7, 'S' : 4, 'T' : 5, 'W' : 11, 'Y' : 7, 'V' : 4, '*' : 1, 'X' : -1, 'B' : 4, 'J' : 3, 'Z' : 4}


def MakeLinksFile(conn_file):
    prot_con_links = {}
    con_prot_links = {}
    with open(conn_file) as file1:
        for line in file1:
            dat = line.strip().split('\t')
            if dat[0] not in con_prot_links:
                con_prot_links[dat[0]] = 0
            con_prot_links[dat[0]] += 1
            prot_con_links[dat[1]] = dat[0]
    return(prot_con_links, con_prot_links)


def MakeBinsDict(conn_file):
    bins_dict = {}
    prot_con_links = {}
    con_prot_links = {}
    with open(conn_file) as file1:
        for line in file1:
            dat = line.strip().split('\t')
            if dat[0] not in bins_dict:
                bins_dict[dat[0]] = {}
            if dat[1] not in bins_dict[dat[0]]:
                bins_dict[dat[0]][dat[1]] = {}
            bins_dict[dat[0]][dat[1]][dat[2]] = ''
            if dat[1] not in con_prot_links:
                con_prot_links[dat[1]] = 0
            con_prot_links[dat[1]] += 1
            prot_con_links[dat[2]] = dat[1]
    return(bins_dict, prot_con_links, con_prot_links)

def MakeBinsConnections_concatenate(prefix, bins):
    out_fasta = open('temporal.fasta','w')
    out_conn = open(prefix + '_bins_to_contigs_connections.txt','w')
    list_bin_files = os.listdir(bins)
    bins_dict = {}
    for bin_file in list_bin_files:
        bins_dict[bin_file] = {}
        with open(bins + bin_file) as file1:
            for line in file1:
                if line.startswith('>'):
                    contig_name = line.rstrip().lstrip('>')
                    bins_dict[bin_file][contig_name] = {}
                    out_conn.write(bin_file + '\t' + contig_name + '\n')
                out_fasta.write(line)
    out_fasta.close()
    out_conn.close()
    return(bins_dict)

def MakeBinsConnections(prefix, bins):
    out_conn = open(prefix + '_bins_to_contigs_connections.txt','w')
    list_bin_files = os.listdir(bins)
    bins_dict = {}
    for bin_file in list_bin_files:
        bins_dict[bin_file] = {}
        with open(bins + bin_file) as file1:
            for line in file1:
                if line.startswith('>'):
                    contig_name = line.rstrip().lstrip('>')
                    bins_dict[bin_file][contig_name] = {}
                    out_conn.write(bin_file + '\t' + contig_name + '\n')
    out_conn.close()
    return(bins_dict)
    
def InfoFromProt(protein_fasta, prot_con_links):
    # Calculate artificial bitscore for each predicted proteins
    prot_weight = {}
    with open(protein_fasta) as file1:
        for line in file1:
            if line.startswith('>'):
                gene = re.search('^>(\S+)\s', line).group(1)
                contig = prot_con_links[gene]
                if contig not in prot_weight:
                    prot_weight[contig] = {}
                prot_weight[contig][gene] = 0
            elif not line.startswith('>'):
                for let in line.strip():
                    if let in sub:
                        prot_weight[contig][gene] += sub[let]
                    else:   
                        print(let)
    return(prot_weight)







