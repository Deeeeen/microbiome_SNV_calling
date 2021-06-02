#!/usr/bin/env python3
import pandas as pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

if __name__=="__main__":
    #### sys.argv[1] file should have format: samples num_sequence
	data = pd.read_csv(sys.argv[1], names=["samples", "num_sequence"], delim_whitespace=True)
	outfile = sys.argv[2]

	## Number of mapped reads bar plot
	plt.figure(1, figsize=(8, 8))
	ind = np.arange(42) 
	width = 0.8
	plt.bar([i+width for i in ind], data["num_sequence"], width)
	plt.xlabel("Sample")
	plt.ylabel("Sequence Counts")
	plt.xticks([i+width for i in ind], data["samples"], rotation=70)
	plt.savefig(outfile + "/num_sequence_bar.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	ind = np.arange(42) 
	width = 0.8
	plt.bar([i+width for i in ind], np.log10(data["num_sequence"]), width)
	plt.xlabel("Sample")
	plt.ylabel("log10(Sequence Counts)")
	plt.xticks([i+width for i in ind], data["samples"], rotation=70)
	plt.savefig(outfile + "/num_sequence_bar_log.png")
	plt.close()

	## Number of mapped reads box plot
	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot(data["num_sequence"])
	plt.ylabel("Sequence Counts")
	plt.xlabel("BH1206")
	plt.savefig(outfile + "/num_sequence_box.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot(np.log10(data["num_sequence"]))
	plt.ylabel("log10(Sequence Counts)")
	plt.xlabel("BH1206")
	plt.savefig(outfile + "/num_sequence_box_log.png")
	plt.close()