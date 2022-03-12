# have to import both, but the second one just gives the extract() to get the sequence at a location.
# if you don't have Bio installed, just have to do "conda install biopython" in terminal
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
import matplotlib.pyplot as plt

# get a record of the file, where the first item from parse() is the largest chromosome
file = 'E. coli Nissle Genome.gbff' 
record = list(SeqIO.parse(file, 'gb'))[0]
record

# if you need the sequence of the whole record
seq = record.seq
len(seq) 

# get a list of all the features from the record, note the first feature is always the whole chromosome
# and the source
features = record.features
features[:10]

# taking the second feature and getting its location
feature = features[1]
location = feature.location
location 

# getting the sequence at that location using extract()
location.extract(record)

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
    
frac_A = num_A/len(seq) * 100
frac_T = num_T/len(seq) * 100
frac_G = num_G/len(seq) * 100
frac_C = num_C/len(seq) * 100

print(num_A)
print(num_T)
print(num_G)
print(num_C)

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
plt.show()
plt.savefig('Distribution of Nucleotide Bases.png')