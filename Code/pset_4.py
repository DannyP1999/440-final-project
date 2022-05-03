# have to import both, but the second one just gives the extract() to get the sequence at a location.
# if you don't have Bio installed, just have to do "conda install biopython" in terminal
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

# get a record of the file, where the first item from parse() is the largest chromosome
file = '../Genome Data/E. coli Nissle Genome.gbff'
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
    
    return genes_df

# load excel spreadsheets with genes of interest from E. coli Nissle (EcN)
# cluster spreadsheets and clusters correspond to more specific families of genes
genes_df_all = pd.read_excel('../Genome Data/Metal Transporter Genes in E. coli Nissle (440).xlsx')
genes_df_cus = pd.read_excel('../Genome Data/cus cluster.xlsx') # copper/silver
genes_df_fec = pd.read_excel('../Genome Data/fec cluster.xlsx') # iron
genes_df_zn = pd.read_excel('../Genome Data/zn cluster.xlsx') # zinc
genes_df_sit = pd.read_excel('../Genome Data/sit cluster.xlsx') # iron/manganese

genes_df_all = gene_df_generator(record, genes_df_all) # all genes listed below
genes_df_cus = gene_df_generator(record, genes_df_cus) # copper/silver transporters
genes_df_fec = gene_df_generator(record, genes_df_fec) # iron transporters
genes_df_zn = gene_df_generator(record, genes_df_zn) # zinc transporters
genes_df_sit = gene_df_generator(record, genes_df_sit) # iron/manganese transporters

## display gene dataframes to check
# display(genes_df_all)
# display(genes_df_cus)
# display(genes_df_fec)
# display(genes_df_zn)
# display(genes_df_sit)

# export gene information to excel
genes_df_all.to_excel('../Genome Data/Genes of Interest Nissle.xlsx')

# ____________________________________________________________________________________________________________
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

    msa_d = msa.digitize(alphabet)

    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa_d, background)
    
    return hmm

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
            return print('Must specify nt or aa sequence.')

        # stop, HMMER time! can't touch this
        sequence1 = alignments[0][1]
        sequence2 = alignments[0][0]
        hmm = stop_HMMER_time([sequence1, sequence2])
        hmm_consensus = hmm.consensus

        alignment_scores.append(alignments[0][2])
        alignment_scores_end.append(alignments[0][4])
        alignment_scores_similarity.append(alignments[0][2]/alignments[0][4])
        consensus_seq.append(hmm_consensus)

    # creating a dictionary and dataframe containing all gene pairs and relevant alignment scores
    gene_pair_names = []
    for i in range(len(gene_pairs)):
        pair_name_OI = genes_df.Gene[gene_pairs[i][0]] + ' + ' + genes_df.Gene[gene_pairs[i][1]]
        gene_pair_names.append(pair_name_OI)

    gene_pair_dict = {'Genes': gene_pair_names, 'Alignment Score': alignment_scores, 
                      'End': alignment_scores_end, 'Similarity': alignment_scores_similarity, 
                      'Consensus': consensus_seq}

    gene_pair_df = pd.DataFrame(data=gene_pair_dict)
    # gene_pair_df
    
    return gene_pair_df


# ____________________________________________________________________________________________________________
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

genes_df_all_paired = gene_pair_df_generator(genes_df_all, nt_or_aa) # all genese listed below
genes_df_cus_paired = gene_pair_df_generator(genes_df_cus, nt_or_aa) # copper/silver transporters
genes_df_fec_paired = gene_pair_df_generator(genes_df_fec, nt_or_aa) # iron transporters
genes_df_zn_paired = gene_pair_df_generator(genes_df_zn, nt_or_aa) # zinc transporters
genes_df_sit_paired = gene_pair_df_generator(genes_df_sit, nt_or_aa) # iron/manganese transporters


# ____________________________________________________________________________________________________________
### Applying HMM
## Align all EcN genes of interest, then build an HMM from the multiple sequence alignment, after which we 
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

genes_df_all = pad_sequences(genes_df_all) # all genes
genes_df_cus = pad_sequences(genes_df_cus) # copper/silver transporters
genes_df_fec = pad_sequences(genes_df_fec) # iron transporters
genes_df_zn = pad_sequences(genes_df_zn) # zinc transporters
genes_df_sit = pad_sequences(genes_df_sit) # iron/manganese transporters

# creating HMM
def create_HMM(genes_df):
    x = []
    for i in range(len(genes_df.Gene)):
        x.append(SeqRecord(Seq(genes_df.Padded.iloc[i]), id=("%s" %genes_df.Gene.iloc[i])))
    align = MultipleSeqAlignment(x)

    align_seq = []
    for i in range(len(align)):
        align_seq.append(str(align[i].seq))

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
        return 0
    
    else:
        ali = hits[0].domains[0].alignment

        target_protein_seq = ali.target_sequence
        target_protein_name = ali.target_name.decode()

        return [target_protein_name, target_protein_seq]

        # Salmonella
faa_salm = "../Genome Data/Salmonella.faa"
# print("\n", faa_salm, "\n")

hmm_all_salm = create_HMM(genes_df_all)
protein_all_salm = apply_HMM(hmm_all_salm, faa_salm)

hmm_cus_salm = create_HMM(genes_df_cus)
protein_cus_salm = apply_HMM(hmm_cus_salm, faa_salm)

hmm_fec_salm = create_HMM(genes_df_fec)
protein_fec_salm = apply_HMM(hmm_fec_salm, faa_salm)

hmm_zn_salm = create_HMM(genes_df_zn)
protein_zn_salm = apply_HMM(hmm_zn_salm, faa_salm)

hmm_sit_salm = create_HMM(genes_df_sit)
protein_sit_salm = apply_HMM(hmm_sit_salm, faa_salm)

salm_proteins_name = [protein_all_salm[0], protein_cus_salm[0],
                      protein_fec_salm[0], protein_zn_salm[0], protein_sit_salm[0]]

salm_proteins_seq = [protein_all_salm[1], protein_cus_salm[1],
                     protein_fec_salm[1], protein_zn_salm[1], protein_sit_salm[1]]

protein_dict_salm = {'Protein Name': salm_proteins_name, 
                     'Protein Sequence': salm_proteins_seq}

protein_df_salm = pd.DataFrame(data=protein_dict_salm)


# Bacillus subtilis
faa_bac_subt = "../Genome Data/Bacillus subtilis.faa"
# print("\n", faa_bac_subt, "\n")

hmm_all_bac_subt = create_HMM(genes_df_all)
protein_all_bac_subt = apply_HMM(hmm_all_bac_subt, faa_bac_subt)

hmm_cus_bac_subt = create_HMM(genes_df_cus)
protein_cus_bac_subt = apply_HMM(hmm_cus_bac_subt, faa_bac_subt)

hmm_fec_bac_subt = create_HMM(genes_df_fec)
protein_fec_bac_subt = apply_HMM(hmm_fec_bac_subt, faa_bac_subt)

hmm_zn_bac_subt = create_HMM(genes_df_zn)
protein_zn_bac_subt = apply_HMM(hmm_zn_bac_subt, faa_bac_subt)

hmm_sit_bac_subt = create_HMM(genes_df_sit)
protein_sit_bac_subt = apply_HMM(hmm_sit_bac_subt, faa_bac_subt)

bac_subt_proteins_name = [protein_all_bac_subt[0], protein_cus_bac_subt[0],
                          protein_fec_bac_subt[0], protein_zn_bac_subt[0], protein_sit_bac_subt[0]]

bac_subt_proteins_seq = [protein_all_bac_subt[1], protein_cus_bac_subt[1],
                         protein_fec_bac_subt[1], protein_zn_bac_subt[1], protein_sit_bac_subt[1]]



protein_dict_bac_subt = {'Protein Name': bac_subt_proteins_name, 
                         'Protein Sequence': bac_subt_proteins_seq}

protein_df_bac_subt = pd.DataFrame(data=protein_dict_bac_subt)


print("\n", faa_salm, "\n")
# display(protein_df_salm)

print("\n", faa_bac_subt, "\n")
# display(protein_df_bac_subt)


# ____________________________________________________________________________________________________________
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
        
    protein_dict_salm_rand = {'Protein_Name': salm_protein_names, 
                              'Protein_Sequence': salm_protein_seqs}

    protein_df_salm_rand = pd.DataFrame(data=protein_dict_salm_rand)


    protein_dict_bac_subt_rand = {'Protein_Name': bac_subt_protein_names, 
                                  'Protein_Sequence': bac_subt_protein_seqs}

    protein_df_bac_subt_rand = pd.DataFrame(data=protein_dict_bac_subt_rand)
    
    return [protein_df_salm_rand, protein_df_bac_subt_rand]


def counts_to_data(counts):
        x = []
        height =[]
        for i in counts:
            x.append(i[0])
            height.append(i[1]/samples)

        return [x, height]

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

    plt.figure(1, figsize = (10, 5))
    bac_subt_labels = ['srfP', 'no hits', 'yscB', 'mntB', 'sppO',
                       'rnY', 'fecD', 'fecE', 'sxzZ', 'yttA', 
                       'ykyA', 'feuV', 'fhuC', 'ykrA', 'pksN',
                       'atpF', 'fecC', 'pspA', 'yoaW', 'spoIIP',
                       'znuC', 'ylqG', 'skiX', 'floT']
    
    # plt.bar(bac_subt_labels[0:(len(counts_bac_subt[0]))], counts_bac_subt[1], color='#FF6103')
    plt.bar(counts_bac_subt[0], counts_bac_subt[1], color='#FF6103')
    plt.xticks(rotation = 45)
    plt.ylabel('Percentage of Bootstrap Replicates')
    plt.legend(['$\it{Bacillus}$ $\it{subtilis}$'])
    fig = plt.xlabel('Gene Identified by HMM')
    
    return fig

samples = 1000
[protein_df_salm_rand, protein_df_bac_subt_rand] = bootstrap(samples, genes_df_all)
# display(protein_df_salm_rand)

fig = plot_bootstrap(protein_df_salm_rand, protein_df_bac_subt_rand)

plt.savefig('../Figures/Salmonella and Bacillus Subtilis Bootstrap.png')

print('done')