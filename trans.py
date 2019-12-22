

def translation(seq, tr_table):
    if tr_table == '1':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '2':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'M','ACA':'T','AAA':'K','AGA':'*',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'*',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '3':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'T','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'T','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'T','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'T','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'M','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '4':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '5':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'M','ACA':'T','AAA':'K','AGA':'S',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'S',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '6':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'Q','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'Q','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '9':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'N','AGA':'S',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'S',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '10':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'C',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '11':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '12':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'S','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '13':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'M','ACA':'T','AAA':'K','AGA':'G',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'G',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '14':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'Y','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'N','AGA':'S',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'S',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '16':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'L','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '21':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'M','ACA':'T','AAA':'N','AGA':'S',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'S',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '22':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'*','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'L','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '23':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'*','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '24':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'S',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'K',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '25':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'G',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '26':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'*','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'A','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '27':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'Q','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'Q','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'A','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '28':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'Q','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'Q','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'A','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '29':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'Y','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'Y','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'A','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '30':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'E','TGA':'*',
    'TTG':'L','TCG':'S','TAG':'E','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'A','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    elif tr_table == '31':
        codontable = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'E','TGA':'W',
    'TTG':'L','TCG':'S','TAG':'E','TGG':'W',
    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',
    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G'
    }
        tr = []
        for i in range(0,len(seq),3):
            if len(seq[i:i+3]) == 3 and seq[i:i+3] in codontable:
                tr.append(codontable[seq[i:i+3]])
            elif len(seq[i:i+3]) == 3 and seq[i:i+3] not in codontable:
                tr.append('X')
        return(''.join(tr))
    else:
        print('Non-exiting translation table')
        exit()
