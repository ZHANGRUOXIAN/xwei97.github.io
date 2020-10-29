# *K*-mer Analysis

Wei Xin in JiangLab, Oct, 2020

## *k*-mer

> In bioinformatics, *k*-mers are subsequences of length ***k*** contained within a biological sequence. *K*-mers are composed of nucleotides (i.e. A, T, G, and C). [[k-mer - Wikipedia](https://en.wikipedia.org/wiki/K-mer)]

## Application

> - Sequence assembly 
> - The detection of genetic islands associated with pathogenicity (**dinucleotide** bias)
> - Genomics-based taxonomy (e.g. **GC-content**)
> - Identify species in metagenomic samples
> - Improve heterologous gene expression (**codon** bias)
> - Create attenuated vaccines (**codon-pair** manipulation)

The frequency of k-mers can be used as a "signature" of sequence.

## Some homemade tools

Recently, I wrote several custom functions that can be used for basic sequence processing and *k*-mers analysis for based on `python3`.

```shell
(base) [weixin@tc6000 ~]$ cd kmer/script/mytools/
(base) [weixin@tc6000 mytools]$ ls
ATGCcontent  DINULC       KmerSpectra  TRANSLATE
COMPLEMENT   KmerCounter  LENGTH       TRINULC
(base) [weixin@tc6000 mytools]$ 
```

The name of each tool corresponds to a Python script and can be invoked independently by `shell` commands.

###  get_fasta()

Load FASTA data. The following tools (functions) all require `get_fasta()` to read the `.fasta` file first. This function can read the `.fasta` file line by line and save the head and sequence as `key` and `value` in a dictionary. 

```python
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
    print("There are %d sequences.\n" % len(fasta))
    return(fasta)
```

### LENGTH

Count DNA bases, i.e. get sequence length.

#### Parameters

| Arguments  | Description                     |
| ---------- | ------------------------------- |
| -h, --help | show this help message and exit |
| -p, -path  | path of your fasta file         |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./LENGTH -h
usage: LENGTH [-h] --path PATH

Count DNA bases, i.e. get sequence length.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
  
(wx) [weixin@tc6000 mytools]$ ./LENGTH -p /public/home/weixin/kmer/fasta/EcoliK12.fasta 

Getting DNA sequence...
There are 1 sequences.

Counting DNA bases...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
length: 4641652 bases
```

The FASTA file contains one sequence of length 4,641,652.

### COMPLEMENT

Get DNA complementary strand.

#### Parameters

| Arguments     | Description                                               |
| ------------- | --------------------------------------------------------- |
| -h, --help    | show this help message and exit                           |
| -p, -path     | path of your fasta file                                   |
| -r, --reverse | whether to reverse complement (default False)             |
| -w, --width   | the number of bases per line output to FASTA (default 60) |
| -o, --outpath | specify the output file (default is to FASTA directory)   |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./COMPLEMENT -h
usage: COMPLEMENT [-h] --path PATH [--reverse REVERSE] [--width WIDTH] [--outpath OUTPATH]

Get DNA complementary strand.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
  --reverse REVERSE, -r REVERSE
                        whether to reverse complement (default False)
  --width WIDTH, -w WIDTH
                        the number of bases per line output to FASTA (default 60)
  --outpath OUTPATH, -o OUTPATH
                        specify the output file (default is to FASTA directory)

(wx) [weixin@tc6000 mytools]$ ./COMPLEMENT -p /public/home/weixin/kmer/fasta/EcoliK12.fasta -r Ture -w 60 -o /public/home/weixin/kmer/fasta/EcoliK12_reverse_complement.fasta

Getting DNA sequence...
There are 1 sequences.

Getting complementary strand...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
OK :)
The FASTA file of the complementary sequence has been saved in '/public/home/weixin/kmer/fasta/EcoliK12_reverse_complement.fasta'.

(wx) [weixin@tc6000 mytools]$ head -5 /public/home/weixin/kmer/fasta/EcoliK12_reverse_complement.fasta
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
GAAAAATACTTACTAAGGCGTTTTTTATTTGGTGATATTTTTTTCAATATCATGCAGCAA
ACGGTGCAACATTGCCGTGTCTCGTTGCTCTAAAAGCCCCAGGCGTTGTTGTAACCAGTC
GACCAGTTTTATGTCATCTGCCACTGCCAGAGTCGTCAGCAATGTCATGGCTCGTTCGCG
TAAAGCTTGCAGTTGATGTTGGTCTGCCGTTGCATCACTTTTCGCCGGTTGTTGTATTAA
```

### TRANSLATE

#### Parameters

| Arguments     | Description                                                 |
| ------------- | ----------------------------------------------------------- |
| -h, --help    | show this help message and exit                             |
| -p, -path     | path of your fasta file                                     |
| -w, --width   | set the number of AA per line output to FASTA (default 60)  |
| -b, --begin   | set translation start location (generally 1/2/3, default 1) |
| -o, --outpath | specify the output file (default is to FASTA directory)     |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./TRANSLATE -h
usage: TRANSLATE [-h] --path PATH [--width WIDTH] [--begin BEGIN] [--outpath OUTPATH]

Translate DNA into AA sequence.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
  --width WIDTH, -w WIDTH
                        set the number of AA per line output to FASTA (default 60)
  --begin BEGIN, -b BEGIN
                        set translation start location (generally 1/2/3, default 1)
  --outpath OUTPATH, -o OUTPATH
                        specify the output file (default is to FASTA directory)

(wx) [weixin@tc6000 mytools]$ ./TRANSLATE -p /public/home/weixin/kmer/fasta/EcoliK12.fasta -w 60 -b 1 -o /public/home/weixin/kmer/fasta/EcoliK12_translated.fasta

Getting DNA sequence...
There are 1 sequences.

Translating DNA into AA sequence...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
OK :)
The FASTA files of the AA sequence has been saved in '/public/home/weixin/kmer/fasta/EcoliK12_translated.fasta'.

(wx) [weixin@tc6000 mytools]$ head -5 /public/home/weixin/kmer/fasta/EcoliK12_translated.fasta
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
SFSF_LQRAICLCVD_KKSV__QLLNWLPAVSKLKFY_LRSLNTLTNIGIAHRQIKITEY
TTSMKRISTTITTTITITTGNGAG_RVQETQKKART_QCGLFFSTKGNEVTTMRVLKFGG
TSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAE
RIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEK
```

### ATGCcontent

#### Description

Count A|T|G|C (1-mers) of DNA sequence (N - unknown base).

#### Parameter

| Arguments  | Description                     |
| ---------- | ------------------------------- |
| -h, --help | show this help message and exit |
| -p, -path  | path of your fasta file         |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./ATGCcontent -h
usage: ATGCcontent [-h] --path PATH

Count A|T|G|C of DNA sequence (N - unknown base).

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file

(wx) [weixin@tc6000 mytools]$ ./ATGCcontent -p /public/home/weixin/kmer/fasta/EcoliK12.fasta 

Getting DNA sequence...
There are 1 sequences.

Counting A|T|G|C of DNA sequence...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
length: 4641652 bases
A: 1142742 (24.6%)
T: 1141382 (24.6%)
G: 1177437 (25.4%)
C: 1180091 (25.4%)
ATGC: 4641652 (100.0%)
N: 0 (0.0%)
```

In this sequence, A/T/G/C accounts for about a quarter, with slightly higher GC content and no unknown nucleotides.

### DINULC

#### Description

Count 2-mers (dinucleotide).

#### Parameters

| Arguments  | Description                     |
| ---------- | ------------------------------- |
| -h, --help | show this help message and exit |
| -p, -path  | path of your fasta file         |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./DINULC -h
usage: DINULC [-h] --path PATH

Count 2-mers.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
(wx) [weixin@tc6000 mytools]$ ./DINULC -p /public/home/weixin/kmer/fasta/EcoliK12.fasta 

Getting DNA sequence...
There are 1 sequences.

Counting 2mers...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
   kmer   count      freq
1    GC  384102  8.275116
11   CG  346793  7.471329
3    TT  339584  7.316018
10   AA  338006  7.282021
5    CA  325327  7.008864
7    TG  322379  6.945352
6    AT  309950  6.677581
15   CC  271821  5.856127
12   GG  270252  5.822325
4    TC  267395  5.760773
8    GA  267384  5.760536
9    AC  256773  5.531932
14   GT  255699  5.508794
0    AG  238013  5.127766
2    CT  236149  5.087608
13   TA  212024  4.567857

OK :)
```

### TRINULC

#### Description

Count 3-mers (trinucleotide).

#### Parameters

| Arguments  | Description                     |
| ---------- | ------------------------------- |
| -h, --help | show this help message and exit |
| -p, -path  | path of your fasta file         |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./TRINULC -h
usage: TRINULC [-h] --path PATH

Count 3-mers.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
(wx) [weixin@tc6000 mytools]$ ./TRINULC -p /public/home/weixin/kmer/fasta/EcoliK12.fasta 

Getting DNA sequence...
There are 1 sequences.

Counting 3mers...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
   kmer   count      freq
57  CGC  115734  2.493381
56  GCG  114670  2.470458
3   TTT  109862  2.366874
34  AAA  108964  2.347527
40  CAG  104850  2.258895
..  ...     ...       ...
19  GGG   47515  1.023666
27  CTC   42746  0.920923
37  GAG   42503  0.915687
39  TAG   27254  0.587162
54  CTA   26770  0.576735

[64 rows x 3 columns]

OK :)
```

### KmerCounter

#### Description

Count *k*-mers.

#### Parameters

| Arguments        | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| -h, --help       | show this help message and exit                              |
| -p, -path        | path of your fasta file                                      |
| -k, --k          | k-mer's k (default 8)                                        |
| -l, --label      | customize a label used to distinguish the results (the label will become legend in the k-mer spectrum |
| -s, --savetofile | whether to save the count for each k-mer to a text file (default True) |
| -o, --outdir     | specify the output directory (default is to FASTA directory) |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./KmerCounter -h
usage: KmerCounter [-h] --path PATH [--k K] [--label LABEL] [--savetofile SAVETOFILE] [--outdir OUTDIR]

Count k-mers.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
  --k K, -k K           k-mer's k (default 8)
  --label LABEL, -l LABEL
                        customize a label used to distinguish the results (the label will become legend in the k-mer spectrum）
  --savetofile SAVETOFILE, -s SAVETOFILE
                        whether to save the count for each k-mer to a text file (default True)
  --outdir OUTDIR, -o OUTDIR
                        specify the output directory (default is to FASTA directory)
                        
(wx) [weixin@tc6000 mytools]$ ./KmerCounter -p /public/home/weixin/kmer/fasta/EcoliK12.fasta -k 8 -l EcoliK12 -s True -o /public/home/weixin/kmer

Getting DNA sequence...
There are 1 sequences.

Counting 8mers...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
           kmer  count      freq     label
580    CGCTGGCG    778  0.016761  EcoliK12
23909  CGCCAGCG    734  0.015813  EcoliK12
13626  CCAGCGCC    726  0.015641  EcoliK12
8748   CGCCAGCA    688  0.014822  EcoliK12
24071  CCGCCAGC    663  0.014284  EcoliK12
...         ...    ...       ...       ...
63192  TTCCAAGG      1  0.000022  EcoliK12
65279  TTCCTAGA      1  0.000022  EcoliK12
65163  TTCCTAGT      1  0.000022  EcoliK12
65338  TTCTAGAT      1  0.000022  EcoliK12
64734  TTGACTAG      1  0.000022  EcoliK12

[65360 rows x 4 columns]

The results of 8-mer counts has been saved in '/public/home/weixin/kmer/EcoliK12_8mer_seq1.txt'
OK :)

(wx) [weixin@tc6000 mytools]$ head -5 /public/home/weixin/kmer/EcoliK12_8mer_seq1.txt
kmer	count	freq	label
CGCTGGCG	778	0.016761	EcoliK12
CGCCAGCG	734	0.015813	EcoliK12
CCAGCGCC	726	0.015641	EcoliK12
CGCCAGCA	688	0.014822	EcoliK12
```

### KmerSpectra

#### Description

Count *k*-mers.

#### Parameters

| Arguments         | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| -h, --help        | show this help message and exit                              |
| -p, -path         | path of your fasta file                                      |
| -k, --k           | k-mer's k (default 8)                                        |
| -l, --label       | customize a label used to distinguish the results (the label will become legend in the k-mer spectrum |
| -s, --savetofile  | whether to save the count for each k-mer to a text file (default True) |
| -o, --outdir      | specify the output directory (default is to FASTA directory) |
| -pdf, --savetopdf | whether to save the spectrum to pdf (default True)           |
| -pp, --pdfpath    | specify the output directory of k-mer spectrum (default is to current working directory) |
| -t, --title       | specify the title of k-mer spectrum                          |
| -c, --col         | specify the colors of points                                 |

#### Examples

```shell
(wx) [weixin@tc6000 mytools]$ ./KmerSpectra  -h
usage: KmerSpectra [-h] --path PATH [--k K] [--label LABEL] [--savetofile SAVETOFILE] [--outdir OUTDIR]
                   [--savetopdf SAVETOPDF] [--pdfpath PDFPATH] [--title TITLE] [--col COL]

Count k-mers.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH, -p PATH  path of your fasta file
  --k K, -k K           k-mer's k (default 8)
  --label LABEL, -l LABEL
                        customize a label used to distinguish the results (the label will become legend in the k-mer spectrum）
  --savetofile SAVETOFILE, -s SAVETOFILE
                        whether to save the count for each k-mer to a text file (default True)
  --outdir OUTDIR, -o OUTDIR
                        specify the output directory of k-mer counts result (default is to FASTA directory)
  --savetopdf SAVETOPDF, -pdf SAVETOPDF
                        whether to save the spectrum to pdf (default True)
  --pdfpath PDFPATH, -pp PDFPATH
                        specify the output directory of k-mer spectrum (default is to current working directory)
  --title TITLE, -t TITLE
                        specify the title of k-mer spectrum
  --col COL, -c COL     specify the colors of points
 
(wx) [weixin@tc6000 mytools]$ ./KmerSpectra -p /public/home/weixin/kmer/fasta/EcoliK12.fasta -k 8 -l EcoliK12 -s True -o /public/home/weixin/kmer -pdf True -pp /public/home/weixin/kmer/pdf/EcoliK12_8mers_spectra.pdf -t 'EcoliK12 8mers spectra' -c '#1f77b4'

Getting DNA sequence...
There are 1 sequences.

Counting 8mers...
>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655, complete genome
           kmer  count      freq     label
580    CGCTGGCG    778  0.016761  EcoliK12
23909  CGCCAGCG    734  0.015813  EcoliK12
13626  CCAGCGCC    726  0.015641  EcoliK12
8748   CGCCAGCA    688  0.014822  EcoliK12
24071  CCGCCAGC    663  0.014284  EcoliK12
...         ...    ...       ...       ...
63192  TTCCAAGG      1  0.000022  EcoliK12
65279  TTCCTAGA      1  0.000022  EcoliK12
65163  TTCCTAGT      1  0.000022  EcoliK12
65338  TTCTAGAT      1  0.000022  EcoliK12
64734  TTGACTAG      1  0.000022  EcoliK12

[65360 rows x 4 columns]

The results of 8-mer counts has been saved in '/public/home/weixin/kmer/EcoliK12_8mer_seq1.txt'
OK :)

Ploting Kmer spectrum...
OK :)
The Kmer spectra has been saved to '/public/home/weixin/kmer/pdf/EcoliK12_8mers_spectra.pdf'
```

An example 8-mer spectrum for E. coli comparing 8-mers' frequency (i.e. multiplicities) with their number of occurrences.

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y86tub3j30yu0ra77r.jpg" alt="image-20201026223929415" style="zoom: 33%;" />



In addition, the *k*-mer spectrum of two sequences can be drawn in one diagram.

```python
# EcoliK12
fasta_path = '/public/home/weixin/kmer/fasta/EcoliK12.fasta'
fasta = get_fasta(fasta_path)
kmer_ecolik12_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                               k=8, savetofile=True, label="EcoliK12")

# EcoliO157
fasta_path = '/public/home/weixin/kmer/fasta/EcoliO157.fasta'
fasta = get_fasta(fasta_path)
kmer_ecolio157_df = KmerCounter(fasta = fasta, fasta_path = fasta_path, 
                                k=8, savetofile=True, label="EcoliO157")

# two Ecoli
import pandas as pd
kmer_ecoli_df  = pd.concat([kmer_ecolik12_df, kmer_ecolio157_df], axis=0)
KmerSpectra(kmer_ecoli_df, savetopdf=True, 
            pdfpath='/public/home/weixin/kmer/pdf/Ecoli_8mer_plot.pdf', 
            title = '8-mer spectrum of ${E.coli}$')
```

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y885tblj30yw0r60w3.jpg" alt="image-20201026224041208" style="zoom:33%;" />

### Dependence

My conda environment:

```shell
(wx) [weixin@tc6000 mytools]$ conda list
# packages in environment at /public/home/weixin/miniconda2/envs/wx:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       1_gnu    conda-forge
ca-certificates           2020.6.20            hecda079_0    conda-forge
certifi                   2020.6.20        py38h924ce5b_2    conda-forge
cycler                    0.10.0                     py_2    conda-forge
dbus                      1.13.16              hb2f20db_0  
expat                     2.2.9                he1b5a44_2    conda-forge
fontconfig                2.13.1            he4413a7_1000    conda-forge
freetype                  2.10.4               he06d7ca_0    conda-forge
glib                      2.66.1               h92f7085_0  
gst-plugins-base          1.14.0               hbbd80ab_1  
gstreamer                 1.14.0               hb31296c_0  
icu                       58.2              hf484d3e_1000    conda-forge
jpeg                      9d                   h516909a_0    conda-forge
kiwisolver                1.2.0            py38hbf85e49_1    conda-forge
lcms2                     2.11                 hbd6801e_0    conda-forge
ld_impl_linux-64          2.35                 h769bd43_9    conda-forge
libblas                   3.8.0               17_openblas    conda-forge
libcblas                  3.8.0               17_openblas    conda-forge
libffi                    3.3                  he6710b0_2  
libgcc-ng                 9.3.0               h5dbcf3e_17    conda-forge
libgfortran-ng            9.3.0               he4bcb1c_17    conda-forge
libgfortran5              9.3.0               he4bcb1c_17    conda-forge
libgomp                   9.3.0               h5dbcf3e_17    conda-forge
liblapack                 3.8.0               17_openblas    conda-forge
libopenblas               0.3.10          pthreads_h4812303_5    conda-forge
libpng                    1.6.37               hed695b0_2    conda-forge
libstdcxx-ng              9.3.0               h2ae2ef3_17    conda-forge
libtiff                   4.1.0                hc7e4089_6    conda-forge
libuuid                   2.32.1            h14c3975_1000    conda-forge
libwebp-base              1.1.0                h516909a_3    conda-forge
libxcb                    1.13              h14c3975_1002    conda-forge
libxml2                   2.9.10               hb55368b_3  
lz4-c                     1.9.2                he1b5a44_3    conda-forge
matplotlib                3.3.2            py38h32f6830_1    conda-forge
matplotlib-base           3.3.2            py38h4d1ce4f_1    conda-forge
ncurses                   6.2                  he1b5a44_2    conda-forge
numpy                     1.19.2           py38hf89b668_1    conda-forge
olefile                   0.46               pyh9f0ad1d_1    conda-forge
openssl                   1.1.1h               h516909a_0    conda-forge
pandas                    1.1.3                    pypi_0    pypi
pcre                      8.44                 he1b5a44_0    conda-forge
pillow                    8.0.0            py38h9776b28_0    conda-forge
pip                       20.2.4                     py_0    conda-forge
pthread-stubs             0.4               h14c3975_1001    conda-forge
pyparsing                 2.4.7              pyh9f0ad1d_0    conda-forge
pyqt                      5.9.2            py38h05f1152_4  
python                    3.8.3                hcff3b4d_2  
python-dateutil           2.8.1                      py_0    conda-forge
python_abi                3.8                      1_cp38    conda-forge
pytz                      2020.1                   pypi_0    pypi
qt                        5.9.7                h5867ecd_1  
readline                  8.0                  he28a2e2_2    conda-forge
scipy                     1.5.3                    pypi_0    pypi
seaborn                   0.11.0                   pypi_0    pypi
setuptools                49.6.0           py38h924ce5b_2    conda-forge
sip                       4.19.13          py38he6710b0_0  
six                       1.15.0             pyh9f0ad1d_0    conda-forge
sqlite                    3.33.0               h4cf870e_1    conda-forge
tk                        8.6.10               hed695b0_1    conda-forge
tornado                   6.0.4            py38h1e0a361_2    conda-forge
wheel                     0.35.1             pyh9f0ad1d_0    conda-forge
xorg-libxau               1.0.9                h14c3975_0    conda-forge
xorg-libxdmcp             1.1.3                h516909a_0    conda-forge
xz                        5.2.5                h516909a_1    conda-forge
zlib                      1.2.11            h516909a_1010    conda-forge
zstd                      1.4.5                h6597ccf_2    conda-forge
```

## Downstream analysis based on *k*-mer frequency

Using the above *k*-mer related tools, I can process *E. coli* K12, *E. coli* O157, *C. crescentus* NA1000 genomes and *H. sapiens* chr21 sequences respectively. Based on the *k*-mer counts list output from `.KmerCounter`, we can analyze the differences and connections between sequences. I used `R` to do the analysis:

- Load packages

```R
library(magrittr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(factoextra)
library(FactoMineR)
library(ggthemes)
```

### Data prepare

#### Read data

Take the 8-MER analysis as an example：

``` R
#### read data ####
ekdf = read.delim2(file = './tables/EcoliK12_8mer_seq1.txt')
eodf = read.delim2(file = './tables/EcoliO157_8mer_seq1.txt')
cndf = read.delim2(file = './tables/caulobacterNA1000_8mer_seq1.txt')
hsdf = read.delim2(file = './tables/HomeSapiens21_8mer_seq1.txt')
dim(hsdf)
hsdf = hsdf[!grepl("N", hsdf$kmer), ]
dim(hsdf) ## [1] 65560     4
hsdf$freq = hsdf$count/sum(hsdf$count) * 100
```

Approximately 14.3% of the *H.sapiens* chr21 sequence is unknown base (N). Therefore, 8-mers containing "N" was eliminated.

#### Merge data

Take the `kmer` and `count` columns from each of the four tables, and merge them into five columns, the first of which is the k-mers after the union.

```R
#### merge data ####
selcol = c('kmer', 'count')
mergedf1 = merge(ekdf[selcol], eodf[selcol], by = 'kmer', all = TRUE, suffixes = c('.ek', '.eo'))
mergedf2 = merge(mergedf1, cndf[selcol], by = 'kmer', all = TRUE)
colnames(mergedf2)[ncol(mergedf2)] %<>% paste0(., '.cn')
mergedf2[is.na(mergedf2)] = 0
mergedf3 = merge(mergedf2, hsdf[selcol], by = 'kmer', all = TRUE)
colnames(mergedf3)[ncol(mergedf3)] %<>% paste0(., '.hs')
mergedf3[is.na(mergedf3)] = 0
mergedf4 = data.frame(
  kmer = mergedf3$kmer,
  apply(mergedf3[, -1], 2, function(x){x/sum(x)*10^5})
)
dim(mergedf4) ## [1] 65560     5
```

The unique 8-mer of each sequence is ≤ 4<sup>8</sup>, after taking the union, all possible 65,560 8-mers are covered.

### Counts distribution

Add up the counts of 8-mer from large to small, and draw the accumulation curves.

```R
# set colors for plot
mycol = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
#### accumulation curve #####
getwd()
accdf = mergedf3[, -1]
acc = NULL
acc.num = NULL
cs_df0 = NULL
for (i in seq_len(ncol(accdf))) {
  c = accdf[, i] %>% sort(., decreasing = TRUE)
  cs_df0 = cbind(cs_df0, cumsum(c))
}
colnames(cs_df0) = c('E.coli K12', 'E.coli O157', 'C.crescentus NA1000', 'H.sapiens chr21')
cs_df1 = data.frame(
  order = 1:nrow(cs_df0),
  cs_df0
)
cs_df2 = melt(cs_df1, id.vars = c('order'), variable.name = 'label', value.name = 'count')
lp1 = ggplot(cs_df2, aes(x=order, y=count, color = label)) +
  geom_line() + 
  xlab('8-mers number') +
  ylab('') + 
  labs(title = 'Accumulation curves of 8-mer counts') +
  scale_color_manual(values = mycol[1:4]) +
  theme_few() +
  theme(legend.position = c(0.75, 0.50))
print(lp1)

lp2 = ggplot(cs_df2[cs_df2$label != 'H.sapiens.chr21', ], aes(x=order, y=count, color = label)) +
  geom_line() + 
  xlab('8-mers number') +
  ylab('') + 
  labs(title = 'Accumulation curves of 8-mer counts') +
  scale_color_manual(values = mycol[1:4]) +
  theme_few() +
  theme(legend.position = c(0.75, 0.25))
print(lp2)
```

The sequence of *H. Sapiens* chr21 is about 10 times the length of these bacterial genomes.

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y84ohiuj31fi0pegpz.jpg" alt="image-20201027095532665" style="zoom:32%;" />

### Calculate frequency

To make the 8-mers of the four sequences comparable, we can divide each count by the total 8-mer of the sequence to get the frequency and then multiply by 10<sup>5</sup>.

```R
#### boxplot ####
bpdf1 = melt(mergedf3, id.vars = c('kmer'), variable.name = 'label', value.name = 'count')
bpdf1$count = log10(bpdf1$count+1)
range(bpdf1$count)
bp1 = ggboxplot(bpdf1, x = 'label', y = "count", 
                color = 'label', bxp.errorbar = T,
                bxp.errorbar.width = 0.3,
                size = 0.3, outlier.size = 1, outlier.shape = 1,
                palette = mycol[1:4]) +
  scale_x_discrete(labels = c('E.coli K12', 'E.coli O157', 'C.crescentus NA1000', 'H.sapiens chr21')) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = '8-mer counts',
       title = 'Distributions of 8-mer counts')
print(bp1)

# calculate frequency
mergedf4 = data.frame(
  kmer = mergedf3$kmer,
  apply(mergedf3[, -1], 2, function(x){x/sum(x)*10^5})
)
dim(mergedf4)

bpdf2 = melt(mergedf4, id.vars = c('kmer'), variable.name = 'label', value.name = 'count')
bpdf2$count = log10(bpdf2$count+1)
range(bpdf2$count)
bp2 = ggboxplot(bpdf2, x = 'label', y = "count", 
                color = 'label', bxp.errorbar = T,
                bxp.errorbar.width = 0.3,
                size = 0.3, outlier.size = 1, outlier.shape = 1,
                palette = mycol[1:4]) +
  scale_x_discrete(labels = c('E.coli K12', 'E.coli O157', 'C.crescentus NA1000', 'H.sapiens chr21')) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = 'log10(8-mer frequency*10e5 + 1)',
       title = 'Distributions of 8-mer frequency')
print(bp2)
```

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y85h5rsj31cy0ri439.jpg" alt="image-20201027124447713" style="zoom:33%;" />

In the counts list, we can find that the two highest outliers of *H. Sapiens* chr21 are corresponding to subsequences poly-T and poly-A, embodying eukaryotic characteristics.

```R
> hsdf = hsdf[order(hsdf$count, decreasing = T), ]
> hsdf[1:2,]
      kmer count      freq         label
2 TTTTTTTT 55795 0.1391803 HomeSapiens21
3 AAAAAAAA 53466 0.1333706 HomeSapiens21
```

In addition, compared with *E. coli*, there are 92 8-mers with count >1000 in *C.* crescentus .

```R
> table(eodf$count > 1000)

FALSE 
65452 
> table(ekdf$count > 1000)

FALSE 
65360 
> table(cndf$count > 1000)

FALSE  TRUE 
64331    92 
```

 It can be seen that most of these sub-sequences are CG bases.

```R
> cndf = cndf[order(cndf$count, decreasing = T), ]
> cndf[cndf$count > 1000, ]
       kmer count                 freq             label
1  CGGCGGCG  2825  0.06987520412216709 caulobacterNA1000
2  CGCCGCCG  2816  0.06965259285239735 caulobacterNA1000
3  GGCGGCGG  2163 0.053500908501326515 caulobacterNA1000
        ··· ··· omit 86 lines here ··· ···
90 GGCGCGGG  1009 0.024957196799740382 caulobacterNA1000
91 GGCGGCCT  1008 0.024932462214210416 caulobacterNA1000
92 GCCCGCGC  1007 0.024907727628680446 caulobacterNA1000
```

### Sequence similarity

#### Clustering

Based on the normalized frequency, we conducted hierarchical clustering to investigated the similarity between sequences (species).

```R
#### Clustering ####
cdf = scale(t(mergedf4[-1]))
row.names(cdf) = c('E.coliK12', 'E.coliO157', 'C.crescentus NA1000', 'H.sapiens chr21')
res <- hcut(cdf, k = 3, stand = TRUE)
fviz_dend(res, k=4,
          cex = 0.5,
          k_colors = mycol[2:4] %>% rev())
```

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y865udij30oq0o4ab8.jpg" alt="image-20201026225431993" style="zoom:33%;" />

Compared with the *H.sapiens* chr21, the three bacterial genomes are also grouped together in a cluster.

#### Principal Components Analysis (PCA)

```R
pcadf = t(mergedf4[-1])
row.names(pcadf) = c('E.coli K12', 'E.coli O157', 'C.crescentus NA1000', 'H.sapiens chr21')
res.pca <- PCA(pcadf,  graph = FALSE)
summary(res.pca)
fviz_pca_ind(res.pca, col.ind = factor(row.names(pcadf), levels = row.names(pcadf)),
             repel = T, geom = c("point", "text"),
             palette = mycol[1:4]) +
  theme(legend.position = "none") +
  scale_shape_manual(values = rep(19, 4))
```

<img src="https://tva1.sinaimg.cn/large/0081Kckwly1gk3y87qlwdj30r20r0gnz.jpg" alt="image-20201026225351661" style="zoom:33%;" />

As we can see from the figures above, the two strains of *E. coli* are the most similar. 

### Dependence

```R
> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.0.2      jiebaR_0.11      jiebaRD_0.1      wordcloud2_0.2.1 ggthemes_4.2.0  
 [6] FactoMineR_2.3   factoextra_1.0.7 reshape2_1.4.4   ggpubr_0.4.0     ggplot2_3.3.2   
[11] magrittr_1.5    

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.0     purrr_0.3.4          haven_2.3.1          lattice_0.20-41     
 [5] carData_3.0-4        colorspace_1.4-1     vctrs_0.3.4          generics_0.0.2      
 [9] htmltools_0.5.0      yaml_2.2.1           rlang_0.4.7          pillar_1.4.6        
[13] foreign_0.8-80       glue_1.4.2           withr_2.3.0          readxl_1.3.1        
[17] lifecycle_0.2.0      plyr_1.8.6           stringr_1.4.0        munsell_0.5.0       
[21] ggsignif_0.6.0       gtable_0.3.0         cellranger_1.1.0     zip_2.1.1           
[25] htmlwidgets_1.5.2    leaps_3.1            labeling_0.3         rio_0.5.16          
[29] forcats_0.5.0        curl_4.3             broom_0.7.1          Rcpp_1.0.5          
[33] scales_1.1.1         backports_1.1.10     flashClust_1.01-2    jsonlite_1.7.1      
[37] scatterplot3d_0.3-41 abind_1.4-5          farver_2.0.3         digest_0.6.25       
[41] hms_0.5.3            stringi_1.5.3        openxlsx_4.2.2       rstatix_0.6.0       
[45] ggrepel_0.8.2        grid_4.0.2           tools_4.0.2          tibble_3.0.3        
[49] cluster_2.1.0        crayon_1.3.4         car_3.0-10           tidyr_1.1.2         
[53] pkgconfig_2.0.3      ellipsis_0.3.1       MASS_7.3-51.6        data.table_1.13.0   
[57] rstudioapi_0.11      R6_2.4.1             compiler_4.0.2      
```

