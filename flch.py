import subprocess
import sys
import os

def CheckSpace(line, log, quiet):
    if ' ' in line:
        if quiet == False:
            print('WARNING: sequence headers contain spaces, this can create problems in classification as diamond outputs header before the first space only.')
        log.write('WARNING: sequence headers contain spaces, this can create problems in classification as diamond outputs header before the first space only.\n')
        return(True)

def CheckDNA(line, log, quiet):
    if len(set(line)) > 5:
        if quiet == False:
            print('WARNING: Are you sure that provided sequences are DNA?')
        log.write('WARNING: Are you sure that provided sequences are DNA?\n')
        return(True)


def CheckDiamond(diamond, db, log, quiet):
    try:
        subprocess.call([diamond, 'help'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    except OSError:
        if quiet == False:
            print('Can not find diamond. Please check whether it is installed or path to binaries is provided')
        log.write('Can not find diamond. Please check whether it is installed or path to binaries is provided' + '\n')
        log.close()
        sys.exit()

def CheckTabFiles(file_to_check, log, quiet):
    try:
        a = open(file_to_check)
        for line in a:
            dat = line.strip().split('\t')
            if len(dat) < 2:
                if quiet == False:
                    print('File ' + file_to_check + ' is not tab-delimited or empty')
                log.write('File ' + file_to_check + ' is not tab-delimited or empty' + '\n')
                log.close()
                sys.exit()
        a.close()
    except OSError:
        if quiet == False:
            print('Cannot find ' + file_to_check)
        log.write('Cannot find ' + file_to_check + '\n')
        log.close()
        sys.exit(1)


def CheckTaxFiles(path_to_taxonomy_files, log, quiet):
    try:
        a = open(path_to_taxonomy_files + 'names.dmp')
        a.close()
    except OSError:
        if quiet == False:
            print('Cannot find names.dmp taxonomy file in the directory ' + path_to_taxonomy_files)
        log.write('Cannot find names.dmp taxonomy file in the directory ' + path_to_taxonomy_files + '\n')
        log.close()
        sys.exit(1)
    try:
        a = open(path_to_taxonomy_files + 'prot.accession2taxid')
        a.close()
    except OSError:
        if quiet == False:
            print('Cannot find prot.accession2taxid taxonomy file in the directory ' + path_to_taxonomy_files)
        log.write('Cannot find prot.accession2taxid taxonomy file in the directory ' + path_to_taxonomy_files + '\n')
        log.close()
        sys.exit(1)
    try:
        a = open(path_to_taxonomy_files + 'nodes.dmp')
        a.close()
    except OSError:
        if quiet == False:
            print('Cannot find nodes.dmp taxonomy file in the directory ' + path_to_taxonomy_files)
        log.write('Cannot find nodes.dmp taxonomy file in the directory ' + path_to_taxonomy_files + '\n')
        log.close()
        sys.exit(1)


def CheckProdigal(prodigal, log, quiet):
    try:
        subprocess.call([prodigal], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    except OSError:
        if quiet == False:
            print('Can not find prodigal. Please check whether it is installed or path to binaries is provided')
        log.write('Can not find prodigal. Please check whether it is installed or path to binaries is provided' + '\n')
        log.close()
        sys.exit()

def CheckFna(fna, log, quiet):
    try:
        a = open(fna)
        a.close()
    except OSError:
        if quiet == False:
            print('Can not find --fna file. Please check whether it is provided')
        log.write('Can not find --fna file. Please check whether it is provided' + '\n')
        log.close()
        sys.exit()
    contigs_info = {}
    space = False
    dna = False
    with open(fna) as file1:
        contigs_names = []
        contigs_names_set = set()
        loc = []
        contig = ''
        for line in file1:
            if line.startswith('>'):
                if space == False:
                    CheckSpace(line.strip(), log, quiet)
                if loc:
                    contigs_info[contig] = {}
                    contigs_info[contig]['length'] = len(''.join(loc))
                    contigs_info[contig]['orfs'] = []
                    loc = []
                contig = line.rstrip().lstrip('>')
                contigs_names.append(contig)
                contigs_names_set.add(contig)
            else:
                if dna == False:
                    dna = CheckDNA(line.strip(), log, quiet)
                loc.append(line.strip())
        contigs_info[contig] = {}
        contigs_info[contig]['length'] = len(''.join(loc))
        contigs_info[contig]['orfs'] = []
    if len(contigs_names) == 0:
        print('Apparently file with DNA sequences is not in FASTA format. Exit...')
        sys.exit()
    if len(contigs_names) != len(contigs_names_set):
        print('Contigs names in provided fasta file are non-unique. Exit...')
        sys.exit()
    else:
        return(contigs_info)

def CheckFna_bins(bins, log, quiet):
    list_bin_files = os.listdir(bins)
    try:
        a = open(bins + list_bin_files[0])
        a.close()
    except OSError:
        if quiet == False:
            print('Can not find directory with bins. Please check whether it is provided')
        log.write('Can not find directory with bins. Please check whether it is provided' + '\n')
        log.close()
        sys.exit()
    contigs_info = {}
    if len(list_bin_files) == 0:
        if quiet == False:
            print('Please check presence of files with bins in the directory ' + bins)
        log.write('Please check presence of files with bins in the directory ' + bins + '\n')
        log.close()
        sys.exit()
    space = False
    dna = False
    for bin_file in list_bin_files:
        with open(bins + bin_file) as file1:
            contigs_names = []
            contigs_names_set = set()
            loc = []
            contig = ''
            for line in file1:
                if line.startswith('>'):
                    if space == False:
                        space = CheckSpace(line.strip(), log, quiet)
                    if loc:
                        contigs_info[contig] = {}
                        contigs_info[contig]['length'] = len(''.join(loc))
                        contigs_info[contig]['orfs'] = []
                        loc = []
                    contig = line.rstrip().lstrip('>')
                    contigs_names.append(contig)
                    contigs_names_set.add(contig)
                else:
                    if dna == False:
                        dna = CheckDNA(line.strip(), log, quiet)
                    loc.append(line.strip())
            contigs_info[contig] = {}
            contigs_info[contig]['length'] = len(''.join(loc))
            contigs_info[contig]['orfs'] = []
        if len(contigs_names) == 0:
            if quiet == False:
                print('Apparently file ' + bin_file + ' with DNA sequences is not in FASTA format. Exit...')
            log.write('Apparently file ' + bin_file + ' with DNA sequences is not in FASTA format. Exit...\n')
            log.close()
            sys.exit()
        if len(contigs_names) != len(contigs_names_set):
            if quiet == False:
                print('Contigs names in ' + bin_file + ' file are non-unique. Exit...')
            log.write('Contigs names in ' + bin_file + ' file are non-unique. Exit...\n')
            log.close()
            sys.exit()
    if quiet == False:
        print(str(len(list_bin_files)) + ' bins for annotation provided')
    log.write(str(len(list_bin_files)) + ' bins for annotation provided\n')
    return(contigs_info)


def CheckFaa(faa, log, quiet):
    try:
        a = open(faa)
        a.close()
    except OSError:
        if quiet == False:
            print('Can not find --faa file. Please check whether it is provided')
        log.write('Can not find --faa file. Please check whether it is provided\n')
        log.close()
        sys.exit()
    space = False
    with open(faa) as file1:
        fl = 0
        for line in file1:
            if line.startswith('>') and fl == 0:
                fl = 1
                if space == False:
                    CheckSpace(line.strip(), log, quiet)
            elif not line.startswith('>') and fl == 1:
                if len(set(line.strip())) < 5 or len(set(line.strip())) > 21:
                    if quiet == False:
                        print('WARNING: Are you sure that the provided --faa file is protein FASTA?')
                    log.write('WARNING: Are you sure that the provided --faa file is protein FASTA?/n')
                return(True)
            else:
                if quiet == False:
                    print('Apparently file with protein sequences is not in FASTA format. Exit...')
                log.write('Apparently file with protein sequences is not in FASTA format. Exit...\n')
                log.close()
                sys.exit()

def CheckDiamondFile(d, log, quiet):
    try:
        a = open(d)
        a.close()
    except OSError:
        if quiet == False:
            print('Can not find -d file. Please check whether it is provided')
        log.write('Can not find -d file. Please check whether it is provided\n')
        log.close()
        sys.exit()
    with open(d) as file1:
        for line in file1:
            dat = line.strip().split('\t')
            if len(dat) > 3:
                return(True)
            else:
                if quiet == False:
                    print('Apparently alignment file is not tab-delimited and/or not in m8 format. Exit...')
                log.write('Apparently alignment file is not tab-delimited and/or not in m8 format. Exit...\n')
                log.close()
                sys.exit()

def runProdigal(prefix, fasta, prodigal, log, quiet):
    protein_gff = prefix + '_orfs.gff'
    protein_fasta = prefix + '_predicted_prot.faa'
    if quiet == False:
        print('CAT run in ORF prediction mode. ORFs prediction is starting...')
    log.write('CAT run in ORF prediction mode. ORFs prediction is starting...\n')
    try:
        subprocess.check_call([prodigal, '-i', fasta, '-a', protein_fasta, '-o', protein_gff, '-p', 'meta', '-q', '-f', 'gff'])
    except subprocess.CalledProcessError:
        if quiet == False:
            print('Prodigal crashed...')
        log.write('Prodigal crashed...\n')
        log.close()
        sys.exit(2)
    prot_con_links = {}
    con_prot_links = {}
    with open(protein_gff) as file1:
        for line in file1:
            if not line.startswith('#'):
                dat = line.strip().split('\t')
                if dat[0] not in con_prot_links:
                    con_prot_links[dat[0]] = 0
                con_prot_links[dat[0]] += 1
    if quiet == False:
        print('ORFs prediction has finished. File with contigs to proteins links is creating...')
    log.write('ORFs prediction has finished. File with contigs to proteins links is creating...' + '\n')

    with open(prefix + '_contigs_to_protein_names_connections.txt','w') as out_file:
        for item in con_prot_links:
            for i in range(con_prot_links[item]):
                prot_con_links[item + '_' + str(i + 1)] = item
                out_file.write(item + '\t' + item + '_' + str(i + 1) + '\n')
    return(prot_con_links, con_prot_links)

def runProdigal_bat(prefix, prodigal, bins_dict, log, quiet):
    protein_gff = prefix + '_orfs.gff'
    protein_fasta = prefix + '_predicted_prot.faa'
    fasta = 'temporal.fasta'
    if quiet == False:
        print('CAT run in ORF prediction mode. ORFs prediction is starting...')
    log.write('CAT run in ORF prediction mode. ORFs prediction is starting...\n')

    try:
        subprocess.check_call([prodigal, '-i', fasta, '-a', protein_fasta, '-o', protein_gff, '-p', 'meta', '-q', '-f', 'gff'])
    except subprocess.CalledProcessError:
        if quiet == False:
            print('Prodigal crashed...')
        log.write('Prodigal crashed...\n')
        log.close()
        sys.exit(2)
    prot_con_links = {}
    con_prot_links = {}
    with open(protein_gff) as file1:
        for line in file1:
            if not line.startswith('#'):
                dat = line.strip().split('\t')
                if dat[0] not in con_prot_links:
                    con_prot_links[dat[0]] = 0
                con_prot_links[dat[0]] += 1
    if quiet == False:
        print('ORFs prediction has finished. File with bins to contigs to proteins links is creating...')
    log.write('ORFs prediction has finished. File with bins to contigs to proteins links is creating...' + '\n')

    with open(prefix + '_bins_to_contigs_to_protein_names_connections.txt','w') as out_file:
        for item in con_prot_links:
            for i in range(con_prot_links[item]):
                prot_con_links[item + '_' + str(i + 1)] = item
            for file_name in bins_dict:
                if item in bins_dict[file_name]:
                    for i in range(con_prot_links[item]):
                        bins_dict[file_name][item][item + '_' + str(i + 1)] = ''
                        out_file.write(file_name + '\t' + item + '\t' + item + '_' + str(i + 1) + '\n')
    return(prot_con_links, bins_dict)



def runDiamond(diamond, prefix, database, protein_fasta, sens, threads, log, quiet):
    diamond_file = prefix + '_alignment_diamond'
    if quiet == False:
        print('Homology search with diamond is starting...')
    log.write('Homology search with diamond is starting...' + '\n')
    try:
        if sens:
            if threads:
                subprocess.check_call([diamond, 'blastp', '-d', database, '-q', protein_fasta, '-k', '0', '--sensitive', '-p', str(threads), '-o', diamond_file, '--quiet'])
            else:
                subprocess.check_call([diamond, 'blastp', '-d', database, '-q', protein_fasta, '-k', '0', '--sensitive', '-o', diamond_file, '--quiet'])
        else:
            if threads:
                subprocess.check_call([diamond, 'blastp', '-d', database, '-q', protein_fasta, '-k', '0', '-p', str(threads), '-o', diamond_file, '--quiet'])
            else:
                subprocess.check_call([diamond, 'blastp', '-d', database, '-q', protein_fasta, '-k', '0', '-o', diamond_file, '--quiet'])
    except subprocess.CalledProcessError:
        if quiet == False:
            print('Something went wrong with diamond. Please check database location')
        log.write('Something went wrong with diamond. Please check database location\n')
        log.close()
        sys.exit()
    if quiet == False:
        print('Homology search with diamond has finished...')
    log.write('Homology search with diamond has finished...' + '\n')

    







    
    
    
    
    
