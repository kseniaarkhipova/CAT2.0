import trans
import random
import os

def fragments_translation(fasta, fragment, prefix, tr_table):
    file_out = open(prefix + '_fragments_6_frames_translations.faa', 'w')
    file_conn = open(prefix + '_contigs_to_protein_names_connections.txt','w')
    con_prot_links = {}
    def seq_check(x):
        return(''.join([{'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'N': 'N', 'R' : random.choice(['A', 'G']), 'Y' : random.choice(['C', 'T']), 'S' : random.choice(['C', 'G']), 'W' : random.choice(['A', 'T']), 'K' : random.choice(['T', 'G']), 'M' : random.choice(['A', 'C']), 'B' : random.choice(['C', 'G', 'T']), 'D' : random.choice(['A', 'G', 'T']), 'H' : random.choice(['A', 'C', 'T']), 'V' : random.choice(['A', 'G', 'C'])}[B] for B in x]))
    def revcompl(x):
        return(''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1]))
    def handling(sequence, name):
        seq = seq_check(sequence)
        count = 1
        if len(seq) < fragment:
            frame_1 = trans.translation(seq, tr_table)
            frame_2 = trans.translation(seq[1:], tr_table)
            frame_3 = trans.translation(seq[2:], tr_table)
            frag_seq_revcomp = revcompl(seq)
            frame_4 = trans.translation(frag_seq_revcomp, tr_table)
            frame_5 = trans.translation(frag_seq_revcomp[1:], tr_table)
            frame_6 = trans.translation(frag_seq_revcomp[2:], tr_table)
            file_out.write('>' + name + '_fragment_1' + '_frame_1' + '\n' + frame_1 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_2' + '\n' + frame_2 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_3' + '\n' + frame_3 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_4' + '\n' + frame_4 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_5' + '\n' + frame_5 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_6' + '\n' + frame_6 + '\n')
            file_conn.write(name + '\t' + name + '_fragment_1' + '\n')
            con_prot_links[name] = 1
        else:
            con_prot_links[name] = 0
            for i in range(0,len(seq),fragment):
                frag_seq = seq[i:i + fragment]
                if len(frag_seq) >= 100:
                    frame_1 = trans.translation(frag_seq, tr_table)
                    frame_2 = trans.translation(frag_seq[1:], tr_table)
                    frame_3 = trans.translation(frag_seq[2:], tr_table)
                    frag_seq_revcomp = revcompl(frag_seq)
                    frame_4 = trans.translation(frag_seq_revcomp, tr_table)
                    frame_5 = trans.translation(frag_seq_revcomp[1:], tr_table)
                    frame_6 = trans.translation(frag_seq_revcomp[2:], tr_table)
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_1' + '\n' + frame_1 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_2' + '\n' + frame_2 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_3' + '\n' + frame_3 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_4' + '\n' + frame_4 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_5' + '\n' + frame_5 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_6' + '\n' + frame_6 + '\n')
                    file_conn.write(name + '\t' + name + '_fragment_' + str(count) + '\n')
                    con_prot_links[name] += 1
                count += 1
    contigs_names = []
    name = ''
    with open(fasta) as file1:
        for line in file1:
            if line.startswith('>') and not name:
                name = line.lstrip('>').strip()
                contigs_names.append(name)
                seq = []
            elif line.startswith('>') and name:
                handling(''.join(seq), name)
                name = line.lstrip('>').strip()
                contigs_names.append(name)
                seq = []
            elif not line.startswith('>'):
                seq.append(line.strip())
        handling(''.join(seq), name)
    file_out.close()
    file_conn.close()
    return(con_prot_links)

def fragments_translation_bins(fragment, prefix, tr_table, bins):
    file_out = open(prefix + '_fragments_6_frames_translations.faa', 'w')
    file_conn = open(prefix + '_bins_to_contigs_connections.txt','w')
    con_prot_links = {}
    def seq_check(x):
        return(''.join([{'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'N': 'N', 'R' : random.choice(['A', 'G']), 'Y' : random.choice(['C', 'T']), 'S' : random.choice(['C', 'G']), 'W' : random.choice(['A', 'T']), 'K' : random.choice(['T', 'G']), 'M' : random.choice(['A', 'C']), 'B' : random.choice(['C', 'G', 'T']), 'D' : random.choice(['A', 'G', 'T']), 'H' : random.choice(['A', 'C', 'T']), 'V' : random.choice(['A', 'G', 'C'])}[B] for B in x]))
    def revcompl(x):
        return(''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1]))
    def handling(sequence, name, con_prot_links, bin_file, bins_dict):
        seq = seq_check(sequence)

        count = 1
        if len(seq) < fragment:
            frame_1 = trans.translation(seq, tr_table)
            frame_2 = trans.translation(seq[1:], tr_table)
            frame_3 = trans.translation(seq[2:], tr_table)
            frag_seq_revcomp = revcompl(seq)
            frame_4 = trans.translation(frag_seq_revcomp, tr_table)
            frame_5 = trans.translation(frag_seq_revcomp[1:], tr_table)
            frame_6 = trans.translation(frag_seq_revcomp[2:], tr_table)
            file_out.write('>' + name + '_fragment_1' + '_frame_1' + '\n' + frame_1 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_2' + '\n' + frame_2 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_3' + '\n' + frame_3 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_4' + '\n' + frame_4 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_5' + '\n' + frame_5 + '\n')
            file_out.write('>' + name + '_fragment_1' + '_frame_6' + '\n' + frame_6 + '\n')
            file_conn.write(bin_file + '\t' + name + '\t' + name + '_fragment_1' + '\n')
            con_prot_links[name] = 1
            bins_dict[bin_file][name][name + '_fragment_1'] = ''
        else:
            con_prot_links[name] = 0
            for i in range(0,len(seq),fragment):
                frag_seq = seq[i:i + fragment]
                if len(frag_seq) >= 100:
                    frame_1 = trans.translation(frag_seq, tr_table)
                    frame_2 = trans.translation(frag_seq[1:], tr_table)
                    frame_3 = trans.translation(frag_seq[2:], tr_table)
                    frag_seq_revcomp = revcompl(frag_seq)
                    frame_4 = trans.translation(frag_seq_revcomp, tr_table)
                    frame_5 = trans.translation(frag_seq_revcomp[1:], tr_table)
                    frame_6 = trans.translation(frag_seq_revcomp[2:], tr_table)
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_1' + '\n' + frame_1 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_2' + '\n' + frame_2 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_3' + '\n' + frame_3 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_4' + '\n' + frame_4 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_5' + '\n' + frame_5 + '\n')
                    file_out.write('>' + name + '_fragment_' + str(count) + '_frame_6' + '\n' + frame_6 + '\n')
                    file_conn.write(bin_file + '\t' + name + '\t' + name + '_fragment_' + str(count) + '\n')
                    con_prot_links[name] += 1
                    bins_dict[bin_file][name][name + '_fragment_' + str(count)] = ''
                count += 1
        return(con_prot_links, bins_dict)
        
        

    list_bin_files = os.listdir(bins)
    bins_dict = {}
    for bin_file in list_bin_files:
        name = ''
        bins_dict[bin_file] = {}
        with open(bins + bin_file) as file1:
            for line in file1:
                if line.startswith('>') and not name:
                    name = line.lstrip('>').strip()
                    bins_dict[bin_file][name] = {}
                    seq = []
                elif line.startswith('>') and name:
                    con_prot_links, bins_dict = handling(''.join(seq), name, con_prot_links, bin_file, bins_dict)
                    name = line.lstrip('>').strip()
                    bins_dict[bin_file][name] = {}
                    seq = []
                elif not line.startswith('>'):
                    seq.append(line.strip())
            con_prot_links, bins_dict = handling(''.join(seq), name, con_prot_links, bin_file, bins_dict)


    file_out.close()
    file_conn.close()
    return(con_prot_links, bins_dict)
