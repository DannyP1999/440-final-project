# ______________________________________________________________________________
# if you don't have Bio installed, just have to do "conda 
# install "pip install biopython" in terminal

# packages to import
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 

import logomaker

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import scipy as sp
import numpy as np
import random
from matplotlib import pyplot as plt

import pyhmmer
pyhmmer.__version__

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from matplotlib.pyplot import figure

from random import randrange
from collections import Counter

# ______________________________________________________________________________
# get a record of the file, where the first item from parse() is the largest 
# chromosome
file = '../Genome Data/Genome and FASTA Files/E. coli Nissle Genome.gbff'
record = list(SeqIO.parse(file, 'gb'))[0]

# if you need the sequence of the whole record
seq = record.seq

## Importing and creating dataframe for relevant Nissle parameters
# takes a gbff genebank file record and excel spreadsheet 
# containing gene names of interest
def gene_df_generator(record, genes_df):
    # initialiize empty gene name and sequence arrays
    gene_names = []
    gene_seqs = []
    gene_tls = []
    gene_products = []
    
    # loop through all the genes of interest
    for i in range(len(genes_df)):
        
        # loop through every feature (gene) of the genome
        for j in range(len(record.features)):
            feature_OI = record.features[j] # current feature of interest

            try:
                if feature_OI.qualifiers['gene'][0] == genes_df.Gene[i]:
                    gene_OI = feature_OI.qualifiers['gene'][0] # gene name
                    
                    # get position of gene in genome
                    location_OI = feature_OI.location 

                    # getting the sequence at that location using extract()
                    seq_record_OI = location_OI.extract(record)
                    gene_seq_OI = seq_record_OI.seq
                    
                    # AA sequence of gene
                    gene_tl_OI = feature_OI.qualifiers['translation'][0] 
                    # gene function description
                    gene_product_OI = feature_OI.qualifiers['product'] 

                    # append gene information to corresponding category
                    gene_products.append(gene_product_OI[0])
                    gene_names.append(gene_OI)
                    gene_seqs.append(gene_seq_OI)
                    gene_tls.append(gene_tl_OI)

            except KeyError as KE: 
                # workaround for key error that results for genes that do 
                # not have gene products listed
                KE = KE

    genes_df['Gene'] = gene_names # input gene name into df
    genes_df['Description'] = gene_products # input gene descirption into df
    genes_df['Sequence'] = gene_seqs # input gene nt sequence into df
    genes_df['Translation'] = gene_tls # input gene aa sequence into df
    
    return genes_df # dataframe with all gene information

# load excel spreadsheets with genes of interest from E. coli Nissle (EcN)
# cluster spreadsheets and clusters correspond to more specific families of genes
genes_df_all = pd.read_excel('../Genome Data/Genes of Interest/Metal Transporter Genes in E. coli Nissle (440).xlsx')
genes_df_cus = pd.read_excel('../Genome Data/Genes of Interest/cus cluster.xlsx') # copper/silver
genes_df_fec = pd.read_excel('../Genome Data/Genes of Interest/fec cluster.xlsx') # iron
genes_df_zn = pd.read_excel('../Genome Data/Genes of Interest/zn cluster.xlsx') # zinc
genes_df_sit = pd.read_excel('../Genome Data/Genes of Interest/sit cluster.xlsx') # iron/manganese

genes_df_all = gene_df_generator(record, genes_df_all) # all genes listed below
genes_df_cus = gene_df_generator(record, genes_df_cus) # copper/silver transporters
genes_df_fec = gene_df_generator(record, genes_df_fec) # iron transporters
genes_df_zn = gene_df_generator(record, genes_df_zn) # zinc transporters
genes_df_sit = gene_df_generator(record, genes_df_sit) # iron/manganese transporters

# export gene information to excel
genes_df_all.to_excel('../Genome Data/Genes of Interest/Genes of Interest Nissle.xlsx')

# ______________________________________________________________________________
## stop... HMMR time
# creates and returns a hidden markov model for gene sequences of interest
def stop_HMMER_time(gene_sequences_OI): # you can't touch this
    # HMMER
    alphabet = pyhmmer.easel.Alphabet.amino()
    
    gene_sequences = []
    
    # loop through all genes of interest and convert to utf_8 format
    for i in range(len(gene_sequences_OI)):
        name = "seq%s" %i
        seq = pyhmmer.easel.TextSequence(name=name.encode('utf_8'), 
                                         sequence=gene_sequences_OI[i])
        gene_sequences.append(seq)
        
    # generate multiple sequence allignment needed to build hmm    
    msa  = pyhmmer.easel.TextMSA(name="msa".encode('utf_8'), 
                                 sequences=gene_sequences)

    # alphabet with multiple sequence alignment
    msa_d = msa.digitize(alphabet)

    # build model with alphabet from multiple sequence alignment
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa_d, background)
    
    return hmm # hidden markov model

# generates all consensus sequences and alignment scores 
# for every combination of genes of interest
def alignment_df_generator(gene_pairs, genes_df, nt_or_aa):
    # initialize empty alignment score arrays
    alignment_scores = []
    alignment_scores_end = []
    alignment_scores_similarity = []
    consensus_seq = []

    for i in range(len(gene_pairs)):
        
        if nt_or_aa == 'nt': # use genetic nucleotide sequence
            alignments = pairwise2.align.globalxx(genes_df.Sequence[gene_pairs[i][0]], 
                                                  genes_df.Sequence[gene_pairs[i][1]])
        elif nt_or_aa == "aa": # use translation amino acid sequence
            alignments = pairwise2.align.globalxx(genes_df.Translation[gene_pairs[i][0]], 
                                                  genes_df.Translation[gene_pairs[i][1]])
        else:
            print('Must specify nt or aa sequence.')
            # amino acids or nucleotides
            return

        # stop, HMMER time! can't touch this
        sequence1 = alignments[0][1]
        sequence2 = alignments[0][0]
        hmm = stop_HMMER_time([sequence1, sequence2])
        hmm_consensus = hmm.consensus

        alignment_scores.append(alignments[0][2])
        alignment_scores_end.append(alignments[0][4])
        alignment_scores_similarity.append(alignments[0][2]/alignments[0][4])
        consensus_seq.append(hmm_consensus)

    # creating a dictionary and dataframe containing all gene pairs and 
    # relevant alignment scores
    gene_pair_names = []
    for i in range(len(gene_pairs)):
        pair_name_OI = genes_df.Gene[gene_pairs[i][0]] + ' + ' + genes_df.Gene[gene_pairs[i][1]]
        gene_pair_names.append(pair_name_OI)

    gene_pair_dict = {'Genes': gene_pair_names, 'Alignment Score': alignment_scores, 
                      'End': alignment_scores_end, 'Similarity': alignment_scores_similarity, 
                      'Consensus': consensus_seq}

    gene_pair_df = pd.DataFrame(data=gene_pair_dict)
    # gene_pair_df
    
    return gene_pair_df # dataframe with all pairs of genes compared


# ______________________________________________________________________________
## Compute alignment scores between pairs of genes
def gene_pair_df_generator(genes_df, nt_or_aa):
    # numbers corresponding to each gene of interest
    gene_nums = (np.linspace(0, len(genes_df)-1, len(genes_df))).astype(int)
    
    # all possible pairs in List using list comprehension + enumerate()
    gene_pairs = [(a, b) for idx, a in enumerate(gene_nums) for b in gene_nums[idx + 1:]]
    gene_pair_df = alignment_df_generator(gene_pairs, genes_df, nt_or_aa)
    gene_pair_ranked_df = gene_pair_df.sort_values(by=['Similarity'], ascending=False)
    
    return gene_pair_ranked_df

# nt_or_aa = "nt"
nt_or_aa = "aa"

# all genese listed below
genes_df_all_paired = gene_pair_df_generator(genes_df_all, nt_or_aa) 
# copper/silver transporters
genes_df_cus_paired = gene_pair_df_generator(genes_df_cus, nt_or_aa) 
# iron transporters
genes_df_fec_paired = gene_pair_df_generator(genes_df_fec, nt_or_aa) 
# zinc transporters
genes_df_zn_paired = gene_pair_df_generator(genes_df_zn, nt_or_aa) 
# iron/manganese transporters
genes_df_sit_paired = gene_pair_df_generator(genes_df_sit, nt_or_aa) 


# ______________________________________________________________________________
### Applying HMM
## Align all EcN genes of interest, then build an HMM from the multiple 
## sequence alignment, after which we 
## can apply the HMM to a sequence database

## Compute alignment scores between pairs of genes
# pads gene sequences in order to ensure all genese of interest are 
# the same length
def pad_sequences(genes_df):
# extract all gene amino acid sequences
    sequences = []
    for i in range(len(genes_df.Translation)):
        sequences.append(genes_df.Translation[i])

    # added "-" to make all genes equal lengths
    longest_length = max(len(s) for s in sequences)
    padded_sequences = [s.ljust(longest_length, '-') for s in sequences]
    records = (SeqRecord(Seq(s)) for s in padded_sequences)

    genes_df['Padded'] = padded_sequences
    
    return genes_df

# dataframes with padded sequences
genes_df_all = pad_sequences(genes_df_all) # all genes
genes_df_cus = pad_sequences(genes_df_cus) # copper/silver transporters
genes_df_fec = pad_sequences(genes_df_fec) # iron transporters
genes_df_zn = pad_sequences(genes_df_zn) # zinc transporters
genes_df_sit = pad_sequences(genes_df_sit) # iron/manganese transporters

# creating HMM for gene dataframe
def create_HMM(genes_df):
    x = []
    for i in range(len(genes_df.Gene)):
        x.append(SeqRecord(Seq(genes_df.Padded.iloc[i]), id=("%s" %genes_df.Gene.iloc[i])))
    align = MultipleSeqAlignment(x)

    align_seq = []
    for i in range(len(align)):
        align_seq.append(str(align[i].seq))

    # build HMM
    hmm = stop_HMMER_time(align_seq)
    
    return hmm


# applying HMM
# faa is protein FASTA file
def apply_HMM(hmm, faa):
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)

    with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)

    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
    hits = pipeline.search_hmm(query=hmm, sequences=sequences)
    
    if len(hits) == 0:
        # sometimes HMM yields no hits, indicating a threshold for homology
        return 0
    
    else:
        ali = hits[0].domains[0].alignment

        target_protein_seq = ali.target_sequence
        target_protein_name = ali.target_name.decode()

        return [target_protein_name, target_protein_seq]


# Salmonella HMM applying
faa_salm = "../Genome Data/Genome and FASTA Files/Salmonella.faa"

# all gene HMM
hmm_all_salm = create_HMM(genes_df_all)
protein_all_salm = apply_HMM(hmm_all_salm, faa_salm)
# copper/silver HMM
hmm_cus_salm = create_HMM(genes_df_cus)
protein_cus_salm = apply_HMM(hmm_cus_salm, faa_salm)
# iron HMM
hmm_fec_salm = create_HMM(genes_df_fec)
protein_fec_salm = apply_HMM(hmm_fec_salm, faa_salm)
# zinc HMM
hmm_zn_salm = create_HMM(genes_df_zn)
protein_zn_salm = apply_HMM(hmm_zn_salm, faa_salm)
# iron/manganese HMM
hmm_sit_salm = create_HMM(genes_df_sit)
protein_sit_salm = apply_HMM(hmm_sit_salm, faa_salm)
# names of HMM-identified proteins
salm_proteins_name = [protein_all_salm[0], protein_cus_salm[0],
                      protein_fec_salm[0], protein_zn_salm[0], protein_sit_salm[0]]
# sequences of HMM-identified protein sequences
salm_proteins_seq = [protein_all_salm[1], protein_cus_salm[1],
                     protein_fec_salm[1], protein_zn_salm[1], protein_sit_salm[1]]
# create dictionary with salmonella information
protein_dict_salm = {'Protein Name': salm_proteins_name, 
                     'Protein Sequence': salm_proteins_seq}
# dataframe with salmonella information
protein_df_salm = pd.DataFrame(data=protein_dict_salm)


# Bacillus subtilis
faa_bac_subt = "../Genome Data/Genome and FASTA Files/Bacillus subtilis.faa"
# all gene HMM
hmm_all_bac_subt = create_HMM(genes_df_all)
protein_all_bac_subt = apply_HMM(hmm_all_bac_subt, faa_bac_subt)
# copper/silver HMM
hmm_cus_bac_subt = create_HMM(genes_df_cus)
protein_cus_bac_subt = apply_HMM(hmm_cus_bac_subt, faa_bac_subt)
# iron HMM
hmm_fec_bac_subt = create_HMM(genes_df_fec)
protein_fec_bac_subt = apply_HMM(hmm_fec_bac_subt, faa_bac_subt)
# zinc HMM
hmm_zn_bac_subt = create_HMM(genes_df_zn)
protein_zn_bac_subt = apply_HMM(hmm_zn_bac_subt, faa_bac_subt)
# iron/manganese HMM
hmm_sit_bac_subt = create_HMM(genes_df_sit)
protein_sit_bac_subt = apply_HMM(hmm_sit_bac_subt, faa_bac_subt)
# names of HMM-identified proteins
bac_subt_proteins_name = [protein_all_bac_subt[0], protein_cus_bac_subt[0],
                          protein_fec_bac_subt[0], protein_zn_bac_subt[0], protein_sit_bac_subt[0]]
# sequences of HMM-identified protein sequences
bac_subt_proteins_seq = [protein_all_bac_subt[1], protein_cus_bac_subt[1],
                         protein_fec_bac_subt[1], protein_zn_bac_subt[1], protein_sit_bac_subt[1]]
# dictionary with bacillus subtilis information
protein_dict_bac_subt = {'Protein Name': bac_subt_proteins_name, 
                         'Protein Sequence': bac_subt_proteins_seq}
# dataframe with bacillus subtilis information
protein_df_bac_subt = pd.DataFrame(data=protein_dict_bac_subt)

# export identified genes to excel in the Figures Folder
protein_df_salm.to_excel('../Figures/Identified Salmonella Transporters.xlsx')
protein_df_salm.to_excel('../Figures/Identified Bacillus Subtillis Transporters.xlsx')


# ______________________________________________________________________________
### BOOTSTRAPPING
# function for generating randomly sampled df for use with bootstrapping
def rand_gene_df_generator(genes_df):
    for i in range(len(genes_df)):
        rand_int = randrange(len(genes_df)) # generate random integer to sample
        if i == 0:
            df = genes_df.iloc[[rand_int]]
            concat_df = df
        else:
            concat_df = pd.concat([concat_df,genes_df.iloc[[rand_int]]])
            
    return concat_df


# bootstrap function
def bootstrap(samples, genes_df):
    
    salm_protein_names = []
    salm_protein_seqs = []
    bac_subt_protein_names = []
    bac_subt_protein_seqs = []
    
    for i in range(samples):
        # print('i = ', i)
        rand_gene_df = rand_gene_df_generator(genes_df)
        rand_hmm = create_HMM(rand_gene_df)

        protein_all_salm = apply_HMM(rand_hmm, faa_salm)
        protein_all_bac_subt = apply_HMM(rand_hmm, faa_bac_subt)

        # if there are no hits based on the randomly sampled df
        # return 'no hits'
        if (protein_all_salm == 0) and (protein_all_bac_subt == 0):
            salm_protein_names.append('no hits')
            salm_protein_seqs.append('no hits')
            bac_subt_protein_names.append('no hits')
            bac_subt_protein_seqs.append('no hits')

        # if only the Salmonella model has no hits
        elif protein_all_salm == 0:
            salm_protein_names.append('no hits')
            salm_protein_seqs.append('no hits')
            bac_subt_protein_names.append(protein_all_bac_subt[0])
            bac_subt_protein_seqs.append(protein_all_bac_subt[1])

        # if only the Bacillus subtilis model has no hits
        elif protein_all_bac_subt == 0:
            bac_subt_protein_names.append('no hits')
            bac_subt_protein_seqs.append('no hits')
            salm_protein_names.append(protein_all_salm[0])
            salm_protein_seqs.append(protein_all_salm[1])

        # if both organisms have hits
        else:
            salm_protein_names.append(protein_all_salm[0])
            salm_protein_seqs.append(protein_all_salm[1])

            bac_subt_protein_names.append(protein_all_bac_subt[0])
            bac_subt_protein_seqs.append(protein_all_bac_subt[1])
    # salmonella randomly selected protein results dictionary
    protein_dict_salm_rand = {'Protein_Name': salm_protein_names, 
                              'Protein_Sequence': salm_protein_seqs}
    # salmonella randomly selected protein results dataframe
    protein_df_salm_rand = pd.DataFrame(data=protein_dict_salm_rand)
    # bacillus subtillis randomly selected protein results
    protein_dict_bac_subt_rand = {'Protein_Name': bac_subt_protein_names, 
                                  'Protein_Sequence': bac_subt_protein_seqs}

    protein_df_bac_subt_rand = pd.DataFrame(data=protein_dict_bac_subt_rand)
    
    return [protein_df_salm_rand, protein_df_bac_subt_rand]

# function for calculating how many times each result is generated from
# each of the thousand bootstrap for generating barplot
def counts_to_data(counts):
        x = []
        height =[]
        for i in counts:
            x.append(i[0])
            height.append(i[1]/samples)

        return [x, height]

# plot bar plot showing bootstrapping results
def plot_bootstrap(protein_df_salm_rand, protein_df_bac_subt_rand):
    z_salm = protein_df_salm_rand.Protein_Name
    z_bac_subt = protein_df_bac_subt_rand.Protein_Name
    counts_salm = Counter(z_salm).most_common()
    counts_bac_subt = Counter(z_bac_subt).most_common()

    counts_salm = counts_to_data(counts_salm)
    counts_bac_subt = counts_to_data(counts_bac_subt)

    plt.figure(0, figsize = (10, 5))
    salm_labels = ['zntB', 'STM0562', 'znuA', 'acrB', 'acrD', 
                   'STM0350', 'STM0351', 'sitA', 'znuB', 'znuC', 
                   'sitC', 'zupT', 'acrF']

    # plt.bar(salm_labels[0:(len(counts_salm[0]))], counts_salm[1])
    plt.bar(counts_salm[0], counts_salm[1])
    plt.title('Ratios of HMM-Identified Genes for 1000 Bootstrap Replicates')
    plt.xticks(rotation = 45)
    plt.ylabel('Percentage of Bootstrap Replicates')
    plt.legend(['$\it{Salmonella}$'])
    plt.xlabel('Gene Identified by HMM')
    plt.savefig('../Figures/Salmonella Bootstrap.png')

    plt.figure(1, figsize = (10, 5))
    # manually determined labels from first bootstrap run used in paper
    bac_subt_labels = ['srfP', 'no hits', 'yscB', 'mntB', 'sppO',
                       'rnY', 'fecD', 'fecE', 'sxzZ', 'yttA', 
                       'ykyA', 'feuV', 'fhuC', 'ykrA', 'pksN',
                       'atpF', 'fecC', 'pspA', 'yoaW', 'spoIIP',
                       'znuC', 'ylqG', 'skiX', 'floT']
    ## manually determined labels from first bootstrap run used in paper, code for using
    ## those labels below
    # plt.bar(bac_subt_labels[0:(len(counts_bac_subt[0]))], counts_bac_subt[1], color='#FF6103')
    plt.bar(counts_bac_subt[0], counts_bac_subt[1], color='#FF6103')
    plt.xticks(rotation = 45)
    plt.title('Ratios of HMM-Identified Genes for 1000 Bootstrap Replicates')
    plt.ylabel('Percentage of Bootstrap Replicates')
    plt.legend(['$\it{Bacillus}$ $\it{subtilis}$'])
    fig = plt.xlabel('Gene Identified by HMM')
    plt.savefig('../Figures/Bacillus Subtilis Bootstrap.png')
    return fig

samples = 1000
[protein_df_salm_rand, protein_df_bac_subt_rand] = bootstrap(samples, genes_df_all)
fig = plot_bootstrap(protein_df_salm_rand, protein_df_bac_subt_rand)



# ______________________________________________________________________________
# initialize number of nucleotide bases 
num_A = 0
num_T = 0
num_G = 0
num_C = 0

for i in range(len(seq)): # loop through every base in the genome sequence
    
    if (seq[i] == 'A'):
        num_A += 1
    
    elif (seq[i] == 'T'):
        num_T += 1
    
    elif (seq[i] == 'G'):
        num_G += 1
    
    elif (seq[i] == 'C'):
        num_C += 1
# calculate percentage of genome for each nucleotide
frac_A = num_A/len(seq) * 100
frac_T = num_T/len(seq) * 100
frac_G = num_G/len(seq) * 100
frac_C = num_C/len(seq) * 100

# creating the dataset
data = {'A':frac_A, 'T':frac_T, 'G':frac_G,'C':frac_C}
bases = list(data.keys())
values = list(data.values())

fig = plt.figure(figsize = (10, 5))

# creating the bar plot
plt.bar(bases, values, color ='maroon', width = 0.4)
plt.xlabel("Nucleotide Base")
plt.ylabel("Percentage of Genome (%)")
plt.title("Distribution of Nucleotide Bases in E. coli Nissle Genome")
plt.savefig('../Figures/Distribution of Nucleotide Bases.png')
plt.show()


# ______________________________________________________________________________
### Onehot Encoding
# Basic structure inspired by: https://machinelearningmastery.com/how-to-one-hot-encode-sequence-data-in-python/
from numpy import argmax

# inputs
alphabet_nt = 'ATCG'
alphabet_aa = 'ACDEFGHIKLMNPQRSTVWY'
# data = 'AACCTGTGAC'
data_nt = genes_df_all.Sequence[0]
data_aa = genes_df_all.Translation[0]

alphabet = alphabet_nt
data = data_nt

# alphabet = alphabet_aa
# data = data_aa

# print('data:', data)

# encoding
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))


# integer encode input data
integer_encoded = [char_to_int[char] for char in data]
# print('data encode:', integer_encoded)


# one hot encode
onehot_encoded = list()
for value in integer_encoded:
    letter = [0 for _ in range(len(alphabet))]
    letter[value] = 1
    onehot_encoded.append(letter)
# print('One Hot encode:', onehot_encoded)


# one hot dataframe
onehot_df = pd.DataFrame([])

onehot_df['Alphabet'] = list(alphabet)

for i in range(len(onehot_encoded)):
    onehot_df[i] = onehot_encoded[i]

# invert encoding
inverted = int_to_char[argmax(onehot_encoded[0])]


# Function-based Dendrogram Generator

def create_square_mat(genes_df, genes_df_paired):
    
    square_mat_df = pd.DataFrame(np.zeros((len(genes_df.Gene),len(genes_df.Gene))),
                                 index = list(genes_df.Gene),
                                 columns = list(genes_df.Gene))

    np.fill_diagonal(square_mat_df.values, 1)

    gene1 = list(np.zeros((len(genes_df_paired.Genes),1)))
    gene2 = list(np.zeros((len(genes_df_paired.Genes),1)))

    for i in range(len(genes_df_paired.Genes)):
        for j in range(len(genes_df.Gene)):
            gene1[i]=genes_df_paired.Genes[i][0:4]
            gene2[i]=genes_df_paired.Genes[i][7:11]
        square_mat_df.loc[gene1[i],gene2[i]] = genes_df_paired.Similarity[i]
        square_mat_df.loc[gene2[i],gene1[i]] = genes_df_paired.Similarity[i]
    
    return square_mat_df

def create_dendogram(square_mat_df, genes_df, gene_OI):
    
    square_mat_df -= 1
    mat = -square_mat_df
    dists = squareform(mat)
    linkage_matrix = linkage(dists, "single")
    figure(figsize=(10,5))
    dend = dendrogram(linkage_matrix, labels=list(genes_df.Gene))
    plt.title("Dendrogram of Genes")
    plt.ylabel("Dissimilarity of Genes")
    plt.xlabel("Genes")
    plt.ylim([0,1])
    plt.savefig('../Figures/Dendrogram of' + gene_OI + '.png')
    plt.show()
    
    return dend


    # All Dendrograms

square_mat_all_df = create_square_mat(genes_df_all, genes_df_all_paired)
dend_all = create_dendogram(square_mat_all_df, genes_df_all, 'all Genes')
# plt.savefig('../Figures/Dendrogram of all Genes.png')

square_mat_cus_df = create_square_mat(genes_df_cus, genes_df_cus_paired)
dend_cus = create_dendogram(square_mat_cus_df, genes_df_cus, 'Copper and Silver Genes')
# plt.savefig('../Figures/Dendrogram of Copper and Silver Genes.png')

square_mat_fec_df = create_square_mat(genes_df_fec, genes_df_fec_paired)
dend_fec = create_dendogram(square_mat_fec_df, genes_df_fec, 'Iron Genes')
# plt.savefig('../Figures/Dendrogram of Iron Genes.png')

square_mat_zn_df = create_square_mat(genes_df_zn, genes_df_zn_paired)
dend_zn = create_dendogram(square_mat_zn_df, genes_df_zn, 'Zinc Genes')
# plt.savefig('../Figures/Dendrogram of Zinc.png')

square_mat_sit_df = create_square_mat(genes_df_sit, genes_df_sit_paired)
dend_sit = create_dendogram(square_mat_sit_df, genes_df_sit, 'Iron and Manganese Genes')
# plt.savefig('../Figures/Dendrogram of Iron and Manganese Genes.png')


print('done')