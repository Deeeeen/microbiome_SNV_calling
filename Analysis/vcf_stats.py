#!/usr/bin/env python3
from optparse import OptionParser
import pandas as pd
import os
import numpy  as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from collections import defaultdict
from matplotlib.collections import LineCollection

def help():
	print("--d1 <filename>                     The input directory")
	print("--d2 <filename>                     The input directory")
	print("--d3 <filename>                     The input directory")
	print("--d4 <filename>                     The input directory")
	print("-o <filename>                       The output directory")
	exit()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename>")
	parser.add_option('--d1', type="string", nargs=1, dest="input1", help="<filename>")
	parser.add_option('--d2', type="string", nargs=1, dest="input2", help="<filename>")
	parser.add_option('--d3', type="string", nargs=1, dest="input3", help="<filename>")
	parser.add_option('--d4', type="string", nargs=1, dest="input4", help="<filename>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if options.help!=None:
		help()
	if options.input1!=None:
		in_dir1 = options.input1
	else:
		raise "No input directory, please include a input directory --d1 <filename>"
	if options.input2!=None:
		in_dir2 = options.input2
	else:
		raise "No input directory, please include a input directory --d2 <filename>"
	if options.input3!=None:
		in_dir3 = options.input3
	else:
		raise "No input directory, please include a input directory --d3 <filename>"
	if options.input4!=None:
		in_dir4 = options.input4
	else:
		raise "No input directory, please include a input directory --d4 <filename>"
	if options.output!=None:
		outfile = options.output
	else:
		outfile = "output"

	data = defaultdict(pd.DataFrame)
	for in_dir in [in_dir1, in_dir2, in_dir3, in_dir4]:
		temp_ref_data = []
		for f in os.listdir(in_dir):
			if f.endswith("_stats"):
				temp_data = pd.DataFrame(data={"sample": [f[:-15]]})
				with open(os.path.join(in_dir, f), "r") as f_in:
					line = f_in.readline()
					while line:
						if line.startswith("SN\t0\tnumber of SNPs:\t"):
							temp_data["num_SNPs"] = int(line.split()[5])
						elif line.startswith("SN\t0\tnumber of MNPs:\t"):
							temp_data["num_MNPs"] = int(line.split()[5])
						elif line.startswith("SN\t0\tnumber of indels:\t"):
							temp_data["num_indels"] = int(line.split()[5])
						elif line.startswith("SN\t0\tnumber of others:\t"):
							temp_data["num_others"] = int(line.split()[5])
						elif line.startswith("SN\t0\tnumber of multiallelic sites:\t"):
							temp_data["num_multi"] = int(line.split()[6])
						elif line.startswith("SN\t0\tnumber of multiallelic SNP sites:\t"):
							temp_data["num_multiSNPs"] = int(line.split()[7])
						line = f_in.readline()
				temp_ref_data.append(temp_data)
		temp_ref_data = pd.concat(temp_ref_data, ignore_index=True)
		temp_ref_data = temp_ref_data.sort_values(['sample']).reset_index(drop=True)
		data[in_dir.split('/')[0]] = temp_ref_data

	## Number of SNPs box plot
	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot([data[k]["num_SNPs"] for k in data.keys()])
	lines = []
	keys = list(data.keys())
	for sam_i in range(42):
		lines.append([[re_i+1, data[keys[re_i]]["num_SNPs"][sam_i]] for re_i in range(len(keys))])
	lc = LineCollection(lines)
	ax.add_collection(lc)
	ax.set_xticklabels(data.keys())
	plt.ylabel("Number of SNVs")
	plt.xlabel("Reference Genome")
	plt.savefig(outfile + "/num_SNPs_box.png")
	plt.close()

	## Number of SNPs bar plot
	plt.figure(1, figsize=(8, 8))
	ind = np.arange(42) 
	width = 0.2
	keys = list(data.keys())
	for re_i in range(len(keys)):
		plt.bar([i+width*re_i for i in ind], data[keys[re_i]]["num_SNPs"], width, label=keys[re_i])
	plt.xlabel("Sample")
	plt.ylabel("Number of SNVs")
	plt.xticks([i+width for i in ind], data["NZ_CP030777.1"]["sample"], rotation=70)
	plt.legend()
	plt.savefig(outfile + "/num_SNPs_bar.png")
	plt.close()


	