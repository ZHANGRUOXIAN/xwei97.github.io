#!/usr/bin/env python
# coding=utf-8

# get DNA sequences
def get_fasta(fasta_path):
    fasta = {}
    print("\nGetting DNA sequence...")
    with open(fasta_path, 'r') as file:
	    for line in file:
	        if line.startswith(">"):
	            name = line.rstrip()
	            fasta[name] = ''
	            continue
	        fasta[name] += line.rstrip()
    print("There are %d sequences." % len(fasta))
    return(fasta)
    
# count DNA bases
def LENGTH(fasta={}, fasta_path=''):
    if fasta == {}: fasta = get_fasta(fasta_path)
    heads = list(fasta.keys())
    seqs = list(fasta.values())
    seqs_num = len(fasta)
    print("\nCounting DNA bases...")
    for i in range(seqs_num):
        print(heads[i])
        print("length: %d bases" % len(seqs[i]))

# count ATGC
def ATGCcontent(fasta={}, fasta_path=''):
    if fasta == {}: fasta = get_fasta(fasta_path)
    heads = list(fasta.keys())
    seqs = list(fasta.values())
    seqs_num = len(fasta)
    print("\nCounting A|T|G|C of DNA sequence...")
    for i in range(seqs_num):
        seq = seqs[i]
        seq_len = len(seq)
        A = seq.count("A")
        T = seq.count("T")
        G = seq.count("G")
        C = seq.count("C")
        ATGC = A+T+G+C
        N = seq.count("N")
        print(heads[i])
        print("length: %d bases" % seq_len)
        print("A: %d (%.1f%%)" % (A, A/seq_len*100))
        print("T: %d (%.1f%%)" % (T, T/seq_len*100))
        print("G: %d (%.1f%%)" % (G, G/seq_len*100))
        print("C: %d (%.1f%%)" % (C, C/seq_len*100))
        print("ATGC: %d (%.1f%%)" % (ATGC, ATGC/seq_len*100))
        print("N: %d (%.1f%%)" % (N, N/seq_len*100))
    
# get complementary strand
def COMPLEMENT(fasta = {}, fasta_path = '', reverse = False, width = 60, outpath = ''):
    if fasta == {}: fasta = get_fasta(fasta_path)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    heads = list(fasta.keys())
    seqs = list(fasta.values())
    seqs_num = len(fasta)
    # default output path
    filename = fasta_path.split("/")[-1]
    if outpath == '': outpath = fasta_path.replace(filename, filename.split(".")[0] + "_complementary.fasta")
    outfile = open(outpath, "w")
    print("\nGetting complementary strand...")
    for i in range(seqs_num):
        print(heads[i])
        seq = seqs[i]
        if reverse: seq = seq[::-1]
        c_seq = ''.join(complement.get(base, base) for base in seq)        
        outfile.write(heads[i]+"\n")
        start=0; end=start+width
        while end < len(c_seq):
            outfile.write(str(c_seq[start:end])+"\n")
            start+=width;end+=width
        while end >= len(c_seq):
            outfile.write(str(c_seq[start:])+"\n")
            break
    outfile.close()
    print("OK :)")
    print("The FASTA file of the complementary sequence has been saved in '%s'." % outpath)

# translate DNA into AA sequence
def TRANSLATE(fasta = {}, fasta_path = '', begin = 1, width = 60, outpath = ''):
    codon_aa = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
    }
    if fasta == {}: fasta = get_fasta(fasta_path)
    heads = list(fasta.keys())
    seqs = list(fasta.values())
    seqs_num = len(fasta)
    # default output path
    filename = fasta_path.split("/")[-1]
    if outpath == '': outpath = fasta_path.replace(filename, filename.split(".")[0] + "_translated.fasta")
    outfile = open(outpath, "w")
    print("\nTranslating DNA into AA sequence...")
    for i in range(seqs_num):
        print(heads[i])
        aa_seq = ''
        seq = seqs[i]
        for j in range(begin-1, len(seq),3):
            if seq[j:j+3] not in codon_aa:
                aa_seq += 'X'  # 'NNN':X, any AA
                break
            if j+3 <= len(seq)-1:
                aa_seq += codon_aa[seq[j:j+3]]
            else: 
                break
        outfile.write(heads[i]+"\n")
        start = 0; end = start + width
        while end < len(aa_seq):
            outfile.write(str(aa_seq[start:end])+"\n")
            start += width; end += width
        while end >= len(aa_seq):
            outfile.write(str(aa_seq[start:])+"\n")
            break
    outfile.close()
    print("OK :)")
    print("The FASTA files of the AA sequence has been saved in '%s'." % outpath)

# count Kmers
def KmerCounter(fasta = {}, fasta_path = '', k = 8, savetofile = False, label = '', outdir = ''):    
    import pandas as pd
    kmers = {}
    if fasta == {}: fasta = get_fasta(fasta_path)
    heads = list(fasta.keys())
    seqs = list(fasta.values())
    seqs_num = len(fasta)
    # default output path
    if savetofile:
        filename = fasta_path.split("/")[-1]
        if outdir == '': 
            outpath = fasta_path.replace(filename, filename.split(".")[0] + '_' + str(k) + "mer")
        else: outpath = outdir + "/" + filename.split(".")[0] + '_' + str(k) + "mer"
    kmer_res = {}
    print("\nCounting %dmers..." % k)
    for i in range(seqs_num):
        seq = seqs[i]
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j+k]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
        df = pd.DataFrame({
            'kmer': list(kmers.keys()), 
            'count': list(kmers.values())
        })
        df['freq'] = df['count']/sum(df['count']) * 100
        if label != '': df['label'] = label
        df_sorted =  df.sort_values(by=['count', 'kmer'], ascending = (False, True))
        if savetofile: 
            outpath = outpath + '_seq'+str(i+1)+'.txt'
            df_sorted.to_csv(outpath, index=False, sep = '\t')
        print(heads[i]) 
        print(df_sorted)
        print()
        kmer_res[heads[i]] = df_sorted
        if savetofile: print("The results of %d-mer counts has been saved in '%s'" % (k, outpath))
    print("OK :)")
    if seqs_num == 1: return(df_sorted)
    else: return(kmer_res)
    
# count 2-mers
def DINULC(fasta = {}, fasta_path = ''):
    KmerCounter(fasta=fasta, fasta_path=fasta_path, k=2, savetofile=False)
    
# count 3-mers
def TRINULC(fasta = {}, fasta_path = ''):
    KmerCounter(fasta=fasta, fasta_path=fasta_path, k=3, savetofile=False)

# Kmer Spectra
def KmerSpectra(kmer_df, savetopdf = True, pdfpath = '', title = "Kmer spectrum", col = ''):
    from collections import Counter
    import matplotlib.pyplot as plt
    import os
    import seaborn as sns
    if col == '':
        pal = sns.color_palette()
        col = pal.as_hex()
    ci = 0
    ax = plt.subplot()
    print("\nPloting Kmer spectrum...")
    for l in set(kmer_df.label):
        sub_df = kmer_df.loc[kmer_df['label'] == l, ]
        count = list(sub_df['count'])
        count_dict = Counter(count)
        x = list(count_dict.keys())
        y = list(count_dict.values())
        ax.scatter(x, y, c = col[ci], s = 15, label = l,
                   alpha = 0.5, edgecolors = 'none')
        ci += 1
    ax.legend()
    ax.grid(alpha = 0.5)
    ax.set_title(title, fontsize=14)  
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Counts')
    print("OK :)")
    if savetopdf:
        if pdfpath == '': 
            pdfpath = os.getcwd() + "/kmer_spectra.pdf"
        ax.figure.savefig(pdfpath)
        print("The Kmer spectra has been saved to " + "'" + pdfpath + "'\n")
    plt.clf()
    plt.close()

# EcoliK12
fasta_path = '/public/home/weixin/kmer/fasta/EcoliK12.fasta'
fasta = get_fasta(fasta_path)
LENGTH(fasta = fasta)
ATGCcontent(fasta = fasta)
COMPLEMENT(fasta = fasta, fasta_path=fasta_path, reverse=False, width=60)
TRANSLATE(fasta = fasta, fasta_path=fasta_path, begin=1, width=60)
DINULC(fasta = fasta)
TRINULC(fasta = fasta)
kmer_ecolik12_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                               k=8, savetofile=True, label="EcoliK12")
KmerSpectra(kmer_ecolik12_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/EcoliK12_8mer_plot.pdf', 
            title = '8-mer spectrum of ${E.coli}$K12')

# EcoliO157
fasta_path = '/public/home/weixin/kmer/fasta/EcoliO157.fasta'
fasta = get_fasta(fasta_path)
LENGTH(fasta = fasta)
ATGCcontent(fasta = fasta)
COMPLEMENT(fasta = fasta, fasta_path=fasta_path, reverse=False, width=60)
TRANSLATE(fasta = fasta, fasta_path=fasta_path, begin=1, width=60)
DINULC(fasta = fasta)
TRINULC(fasta = fasta)
kmer_ecolio157_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                                k=8, savetofile=True, label="EcoliO157")
KmerSpectra(kmer_ecolio157_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/EcoliO157_8mer_plot.pdf', 
            title = '8-mer spectrum of ${E.coli}$O157')

# caulobacterNA1000
fasta_path = '/public/home/weixin/kmer/fasta/caulobacterNA1000.fasta'
fasta = get_fasta(fasta_path)
LENGTH(fasta = fasta)
ATGCcontent(fasta = fasta)
COMPLEMENT(fasta = fasta, fasta_path=fasta_path, reverse=False, width=60)
TRANSLATE(fasta = fasta, fasta_path=fasta_path, begin=1, width=60)
DINULC(fasta = fasta)
TRINULC(fasta = fasta)
kmer_cbna1000_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                               k=8, savetofile=True, label="caulobacterNA1000")
KmerSpectra(kmer_cbna1000_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/caulobacterNA1000_kmerplot.pdf', 
            title = '8-mer spectrum of ${caulobacter}$NA1000')

# two Ecoli
import pandas as pd
kmer_ecolik_df  = pd.concat([kmer_ecolik12_df, kmer_ecolio157_df], axis=0)
KmerSpectra(kmer_ecolik_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/Ecoli_8mer_plot.pdf', 
            title = '8-mer spectrum of ${E.coli}$')

# Human 21
fasta_path = '/public/home/weixin/kmer/fasta/HomeSapiens21.fasta'
fasta = get_fasta(fasta_path)
LENGTH(fasta = fasta)
ATGCcontent(fasta = fasta)
COMPLEMENT(fasta = fasta, fasta_path=fasta_path, reverse=False, width=60)
TRANSLATE(fasta = fasta, fasta_path=fasta_path, begin=1, width=60)
DINULC(fasta = fasta)
TRINULC(fasta = fasta)
kmer_human21_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                              k=8, savetofile=True, label="HomeSapiens21")
KmerSpectra(kmer_human21_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/HomeSapiens21_kmerplot.pdf', 
            title = '8-mer spectrum of ${Home Sapiens}$ chr21')

# 4 sequences
kmer_df = pd.concat([kmer_ecolik12_df, kmer_ecolio157_df, 
                     kmer_cbna1000_df, kmer_human21_df], 
                    axis=0)
KmerSpectra(kmer_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/all_kmerplot.pdf', 
            title = '8-mer spectrum')



