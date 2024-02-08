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

class MultiThread(threading.Thread):

	def __init__(self, r_df, l_df, o) -> None:
		super().__init__()
		self.read_d = r_df
		self.library_d = l_df
		self.out = o


	def run(self):
		for indexx, roww in self.read_d.iterrows():
			if len(roww['Read sequence']) >= 15:
				data = check_output(["bowtie", "-a", "-p", "4", "-v", "3", "example/GS115/GS115", "-c", "%s"%(roww['Read sequence'])], universal_newlines=True)
				if len(data) != 0:
					data = io.StringIO(data)
					dataframe = pd.read_csv(data, sep='\t', header=None, usecols=[0, 1, 2, 3, 4, 5, 6, 7], engine='python')
					for i, r in dataframe.iterrows():
						genomic_location = str(r[2]) + ':' + str(r[3])
						df = self.library_d.loc[(self.library_d['Genomic location'] == genomic_location) & (self.library_d['Strand'] == r[1])]
						if df.shape[0] != 0:
							for p, q in df.iterrows():
								a = int(roww['read_count'])
								if p not in self.out:
									self.out[p] = a
								else:
									self.out[p] += a


outputs = {}
threads = {}

library_df = pd.read_csv('example/Results/sample1_ExactMatched_reads.csv', header=0, sep=',')
read_df = pd.read_csv('example/Results/sample1_unmatched_reads.csv', header=0, sep=',')

read_df = read_df.dropna()
read_df = read_df.reset_index(drop=True)

read_splits = np.array_split(read_df, 64)

for i in range(0,64):
	outputs[i] = {}
	t = MultiThread(read_splits[i], library_df, outputs[i])
	threads[i] = t
	t.start()

for i in range(0,64):
	threads[i].join()

for i in range(0,64):
	for key in range(0, library_df.shape[0] - 1):
		if key in outputs[i]:
			library_df.at[key, 'sample1_sgRNA_count'] += outputs[i][key]

library_df.to_csv('example/Results/Final_sgRNA_count.csv')


