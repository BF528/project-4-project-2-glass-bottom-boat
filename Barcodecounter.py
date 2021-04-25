#!/usr/bin/env python
# coding: utf-8

import gzip
import operator
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from pylab import savefig
from collections import Counter
import matplotlib.pyplot as plt

sample1 = []
sample2 = []
sample3 = []

with gzip.open("/projectnb/bf528/project_4_scrnaseq/fastq/SRR3879604/SRR3879604_1_bc.fastq.gz", "rt") as handle1:
    for record in SeqIO.parse(handle1, "fastq"):
        if len(record.seq) == 25:
            sample1.append(record.seq)

with gzip.open("/projectnb/bf528/project_4_scrnaseq/fastq/SRR3879605/SRR3879605_1_bc.fastq.gz", "rt") as handle2:
    for record in SeqIO.parse(handle2, "fastq"):
        if len(record.seq) == 25:
            sample2.append(record.seq)

with gzip.open("/projectnb/bf528/project_4_scrnaseq/fastq/SRR3879606/SRR3879606_1_bc.fastq.gz", "rt") as handle3:
    for record in SeqIO.parse(handle3, "fastq"):
        if len(record.seq) == 25:
            sample3.append(record.seq)

# from SeqIO format to string
s1_str = []
s2_str = []
s3_str = []

for line in sample1:
    s1_str.append(str(line))
for line in sample2:
    s2_str.append(str(line))
for line in sample3:
    s3_str.append(str(line))

# keep only the barcode, remove UMI
barcodes1 = []
barcodes2 = []
barcodes3 = []

for line in s1_str:
    barcodes1.append(line[:19])
for line in s2_str:
    barcodes2.append(line[:19])
for line in s3_str:
    barcodes3.append(line[:19])

# Counter gives us how many times each barcode is in each file
c1 = Counter(barcodes1)
c2 = Counter(barcodes2)
c3 = Counter(barcodes3)

# write counters to files
with open("SRR3879604_barcode_counts.txt", "w") as f:
    for key in c1.keys():
        f.write("%s\t%s\n"%(key,c1[key]))
with open("SRR3879605_barcode_counts.txt", "w") as f:
    for key in c2.keys():
        f.write("%s\t%s\n"%(key,c2[key]))
with open("SRR3879606_barcode_counts.txt", "w") as f:
    for key in c3.keys():
        f.write("%s\t%s\n"%(key,c3[key]))

# CALCULATE ECDF
# to find the whitelist threshold

# counter to pandas dataframe
df1 = pd.DataFrame.from_dict(c1, orient="index", columns = ['count'])
df2 = pd.DataFrame.from_dict(c2, orient="index", columns = ['count'])
df3 = pd.DataFrame.from_dict(c3, orient="index", columns = ['count'])

# sort values
df1 = df1.sort_values(by='count', ascending=False)
df2 = df2.sort_values(by='count', ascending=False)
df3 = df3.sort_values(by='count', ascending=False)

# plot ecdf
fig1 = sns.ecdfplot(data=df1,x='count')
fig2 = sns.ecdfplot(data=df2,x='count')
fig3 = sns.ecdfplot(data=df3,x='count')

# save figures
figure1 = fig1.get_figure()
figure1.savefig("ecdf_s1.png")
figure2 = fig2.get_figure()
figure2.savefig("ecdf_s2.png")
figure3 = fig3.get_figure()
figure3.savefig("ecdf_s3.png")

