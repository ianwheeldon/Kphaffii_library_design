#!/usr/bin/env python

import threading
import time
start_time = time.time()
from Bio import SeqIO
from collections import Counter
import pandas as pd
import gffutils
import time
import itertools
from subprocess import Popen, PIPE
from subprocess import check_output
import io
from io import StringIO
import numpy as np

fq_dict = SeqIO.index("example/Results/sample1_5ptrimmed_3ptrimmed.fastq", "fastq")
read_list = []
length = []
for i in fq_dict.keys():
       read_list.append(str(fq_dict[i].seq))
       length.append(len(str(fq_dict[i].seq)))

read_list.sort()
counter_dict = Counter(read_list)
print(Counter(length), flush=True)

df = pd.DataFrame.from_dict(counter_dict, orient='index', columns=['read_count'])
df.reset_index(inplace=True)
read_df = df.rename(columns = {'index':'Read sequence'})
read_df = read_df.dropna()
read_df = read_df.reset_index(drop=True)

library_df = pd.read_csv('../Library_design/example/Results/BEST_LIBRARY.csv', header=0, sep=',')

for index, row in library_df.iterrows():
	data = read_df.loc[(read_df['Read sequence'] == row['Target sequence'])]
	if data.shape[0] > 0:
		library_df.at[index, 'sample1_sgRNA_count'] = int(data['read_count'])
		read_df.drop(data.index[0], inplace=True)
	else:
		library_df.at[index, 'sample1_sgRNA_count'] = 0

read_df.to_csv('example/Results/sample1_unmatched_reads.csv')
library_df.to_csv('example/Results/sample1_ExactMatched_reads.csv')

exact_matching_count = list(library_df['sample1_sgRNA_count'])
print('sample1_sgRNA_count', sum(exact_matching_count), flush=True)
