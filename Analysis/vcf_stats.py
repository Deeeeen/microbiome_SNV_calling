#!/usr/bin/env python3
from optparse import OptionParser
import pandas as pd
import os
import numpy  as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.collections import LineCollection

def help():
	print("--d1 <filename>                     The input directory")
	print("--d2 <filename>                     The input directory")
	print("--d3 <filename>                     The input directory")
	print("--d4 <filename>                     The input directory")
	print("-o <filename>                       The output directory")
	exit()

def countSNVs(df1, df2):
	pos_ls = list(set(df1["pos"].tolist()) | set(df2["pos"].tolist()))
	pos_ls.sort()
	num_SNVs = 0
	for pos in pos_ls:
		row1 = df1.loc[df1["pos"] == pos]
		row2 = df2.loc[df2["pos"] == pos]
		if (len(row1) == 0 and len(row2) == 1):
			REF = row2["REF"].tolist()[0]
			ALTs = row2["ALT"].tolist()[0].split(",")
			GT1 = REF
			d2_GT = row2["GT"].tolist()[0]
			if str(d2_GT) == "0":
				GT2 = REF
			else:
				GT2 = ALTs[int(d2_GT)-1]
			if GT1 != GT2 and len(GT1) == len(REF) and len(GT2) == len(REF):
				num_SNVs += 1
		elif (len(row1) == 1 and len(row2) == 0):
			REF = row1["REF"].tolist()[0]
			ALTs = row1["ALT"].tolist()[0].split(",")
			GT2 = REF
			d1_GT = row1["GT"].tolist()[0]
			if str(d1_GT) == "0":
				GT1 = REF
			else:
				GT1 = ALTs[int(d1_GT)-1]
			if GT1 != GT2 and len(GT1) == len(REF) and len(GT2) == len(REF):
				num_SNVs += 1
		elif (len(row1) == 1 and len(row2) == 1):
			REF = row1["REF"].tolist()[0]
			ALTs = row1["ALT"].tolist()[0].split(",")
			d1_GT = row1["GT"].tolist()[0]
			if str(d1_GT) == "0":
				GT1 = REF
			else:
				GT1 = ALTs[int(d1_GT)-1]
			REF = row2["REF"].tolist()[0]
			ALTs = row2["ALT"].tolist()[0].split(",")
			d2_GT = row2["GT"].tolist()[0]
			if str(d2_GT) == "0":
				GT2 = REF
			else:
				GT2 = ALTs[int(d2_GT)-1]
			if GT1 != GT2 and len(GT1) == len(REF) and len(GT2) == len(REF):
				num_SNVs += 1
		else:
			print("\nERROR\n")
			print(len(row1))
			print(row1)
			print(len(row2))
			print(row2)
	return num_SNVs

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
	data_rg = defaultdict(pd.DataFrame)
	T0 = ["B1_A1", "B1_A10", "B1_A11", "B1_A12", "B1_A13", "B1_A14", "B1_A15", 
		  "B1_A16", "B1_A17", "B1_A18", "B1_A19", "B1_A2", "B1_A20", "B1_A21",
		  "B1_A3", "B1_A4", "B1_A5", "B1_A6", "B1_A7", "B1_A8", "B1_A9"]
	T1 = ["B1_B1", "B1_B10", "B1_B11", "B1_B12", "B1_B13", "B1_B14", "B1_B15",
		  "B1_B16", "B1_B17", "B1_B18", "B1_B19", "B1_B2", "B1_B20", "B1_B21",
		  "B1_B3", "B1_B4", "B1_B5", "B1_B6", "B1_B7", "B1_B8", "B1_B9"]

	#### all XXX_GT files should have format: ref pos REF ALT GT
	for in_dir in [in_dir1, in_dir2, in_dir3, in_dir4]:
		samples = []
		num_SNVs = []
		samples0_rg = []
		samples1_rg = []
		num_SNVs0_rg = []
		num_SNVs1_rg = []
		for i in range(len(T0)):
			f1 = T0[i]+"_filtered_GT"
			temp_data1 = pd.read_csv(os.path.join(in_dir, f1), names=["ref", "pos", "REF", "ALT", "GT"], delim_whitespace=True)
			f2 = T1[i]+"_filtered_GT"
			temp_data2 = pd.read_csv(os.path.join(in_dir, f2), names=["ref", "pos", "REF", "ALT", "GT"], delim_whitespace=True)
			num_SNVs.append(countSNVs(temp_data1, temp_data2))
			samples.append(T1[i][4:])
			num_SNVs0_rg.append(len(temp_data1))
			num_SNVs1_rg.append(len(temp_data2))
		data[in_dir.split('/')[0]] = pd.DataFrame(data={"samples": samples,
														"num_SNVs": num_SNVs})
		data_rg[in_dir.split('/')[0]] = pd.DataFrame(data={"samples": T0+T1,
														"num_SNVs": num_SNVs0_rg+num_SNVs1_rg})

	## Number of SNPs box plot
	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot([data[k]["num_SNVs"] for k in data.keys()])
	lines = []
	keys = list(data.keys())
	for sam_i in range(21):
		lines.append([[re_i+1, data[keys[re_i]]["num_SNVs"][sam_i]] for re_i in range(len(keys))])
	lc = LineCollection(lines, alpha=0.5)
	ax.add_collection(lc)
	ax.set_xticklabels(data.keys())
	plt.ylabel("Number of SNVs")
	plt.xlabel("Reference Genome")
	plt.savefig(outfile + "/num_SNPs_box.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot([np.log10(data[k]["num_SNVs"]) for k in data.keys()])
	lines = []
	keys = list(data.keys())
	for sam_i in range(21):
		lines.append([[re_i+1, np.log10(data[keys[re_i]]["num_SNVs"][sam_i])] for re_i in range(len(keys))])
	lc = LineCollection(lines, alpha=0.5)
	ax.add_collection(lc)
	ax.set_xticklabels(data.keys())
	plt.ylabel("log10(Number of SNVs)")
	plt.xlabel("Reference Genome")
	plt.savefig(outfile + "/num_SNPs_box_log.png")
	plt.close()

	## Number of SNPs bar plot
	plt.figure(1, figsize=(8, 8))
	ind = np.arange(21) 
	width = 0.2
	keys = list(data.keys())
	for re_i in range(len(keys)):
		plt.bar([i+width*re_i for i in ind], data[keys[re_i]]["num_SNVs"], width, label=keys[re_i])
	plt.xlabel("Sample")
	plt.ylabel("Number of SNVs")
	plt.xticks([i+width for i in ind], data["NZ_CP030777.1"]["samples"], rotation=70)
	plt.legend()
	plt.savefig(outfile + "/num_SNPs_bar.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	ind = np.arange(21) 
	width = 0.2
	keys = list(data.keys())
	for re_i in range(len(keys)):
		plt.bar([i+width*re_i for i in ind], np.log10(data[keys[re_i]]["num_SNVs"]), width, label=keys[re_i])
	plt.xlabel("Sample")
	plt.ylabel("log10(Number of SNVs)")
	plt.xticks([i+width for i in ind], data["NZ_CP030777.1"]["samples"], rotation=70)
	plt.legend()
	plt.savefig(outfile + "/num_SNPs_bar_log.png")
	plt.close()

	## ALL Number of SNPs box plot
	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot([data_rg[k]["num_SNVs"] for k in data_rg.keys()])
	lines = []
	keys = list(data_rg.keys())
	for sam_i in range(42):
		lines.append([[re_i+1, data_rg[keys[re_i]]["num_SNVs"][sam_i]] for re_i in range(len(keys))])
	lc = LineCollection(lines, alpha=0.5)
	ax.add_collection(lc)
	ax.set_xticklabels(data_rg.keys())
	plt.ylabel("Number of SNVs")
	plt.xlabel("Reference Genome")
	plt.savefig(outfile + "/all_num_SNPs_box.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	fig, ax = plt.subplots()
	ax.boxplot([np.log10(data_rg[k]["num_SNVs"]) for k in data_rg.keys()])
	lines = []
	keys = list(data_rg.keys())
	for sam_i in range(42):
		lines.append([[re_i+1, np.log10(data_rg[keys[re_i]]["num_SNVs"][sam_i])] for re_i in range(len(keys))])
	lc = LineCollection(lines, alpha=0.5)
	ax.add_collection(lc)
	ax.set_xticklabels(data_rg.keys())
	plt.ylabel("log10(Number of SNVs)")
	plt.xlabel("Reference Genome")
	plt.savefig(outfile + "/all_num_SNPs_box_log.png")
	plt.close()

	## ALL Number of SNPs bar plot
	plt.figure(1, figsize=(8, 8))
	ind = np.arange(42) 
	width = 0.2
	keys = list(data_rg.keys())
	for re_i in range(len(keys)):
		plt.bar([i+width*re_i for i in ind], data_rg[keys[re_i]]["num_SNVs"], width, label=keys[re_i])
	plt.xlabel("Sample")
	plt.ylabel("Number of SNVs")
	plt.xticks([i+width for i in ind], data_rg["NZ_CP030777.1"]["samples"], rotation=70)
	plt.legend()
	plt.savefig(outfile + "/all_num_SNPs_bar.png")
	plt.close()

	plt.figure(1, figsize=(8, 8))
	ind = np.arange(42) 
	width = 0.2
	keys = list(data_rg.keys())
	for re_i in range(len(keys)):
		plt.bar([i+width*re_i for i in ind], np.log10(data_rg[keys[re_i]]["num_SNVs"]), width, label=keys[re_i])
	plt.xlabel("Sample")
	plt.ylabel("log10(Number of SNVs)")
	plt.xticks([i+width for i in ind], data_rg["NZ_CP030777.1"]["samples"], rotation=70)
	plt.legend()
	plt.savefig(outfile + "/all_num_SNPs_bar_log.png")
	plt.close()