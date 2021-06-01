#!/usr/bin/env python3
from optparse import OptionParser
import pandas as pd
import os
import numpy  as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

def help():
	print("-d <filename>                     The input directory")
	print("-o <filename>                     The output png file")
	exit()

if __name__=="__main__":
	## read in arguments
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename>")
	parser.add_option('-d', type="string", nargs=1, dest="input", help="<filename>")
	parser.add_option('--minDep', type="int", nargs=1, dest="minDep", help="Minimum Depth")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if options.help!=None:
		help()
	if options.input!=None:
		in_dir = options.input
	else:
		raise "No input directory, please include a input directory -d <filename>"
	if options.minDep!=None:
		minDep = options.minDep
	else:
		raise "Please provide a minDep: --minDep minimum depth"
	if options.output!=None:
		outfile = options.output
	else:
		outfile = "output"

	## read in files
	num_reads = defaultdict()
	mapped_reads = defaultdict()
	depth = defaultdict(pd.DataFrame)
	max_depth = 0
	max_pos = 0
	coverage_list = []
	coverage_minDep_list = []
	acc_coverage = []
	acc_coverage_minDep = []
	total_depth = []
	total_depth_minDep = []
	for f in os.listdir(in_dir):
		if f.endswith("_stats"):
			with open(os.path.join(in_dir, f), "r") as f_in:
				line = f_in.readline()
				while line:
					if line.startswith("SN\tsequences:\t"):
						num_reads[f[:-6]] = int(line.split()[2])
					elif line.startswith("SN\treads mapped:\t"):
						mapped_reads[f[:-6]] = int(line.split()[3])
					line = f_in.readline()
		elif f.endswith("_depth"):
			data = pd.read_csv(os.path.join(in_dir, f), names=["ref", "pos", "depth"], delim_whitespace=True)
			data["pos"] = data["pos"]/1000000
			depth[f[:-6]] = data
			coverage_list.append(float(len(data.loc[data["depth"] > 0])/len(data)))
			coverage_minDep_list.append(float(len(data.loc[data["depth"] >= minDep])/len(data)))
			if max(depth[f[:-6]]["depth"]) > max_depth:
				max_depth = max(depth[f[:-6]]["depth"])
			if max(depth[f[:-6]]["pos"]) > max_pos:
				max_pos = max(depth[f[:-6]]["pos"])
			if len(total_depth) == 0:
				total_depth_minDep = [i if i >= minDep else 0 for i in data["depth"]]
				total_depth = data["depth"]
			else:
				total_depth_minDep = [x + y for x, y in zip(total_depth_minDep, [i if i >= minDep else 0 for i in data["depth"]])]
				total_depth = [x + y for x, y in zip(total_depth, data["depth"])]
			acc_coverage.append(float(np.count_nonzero(total_depth)/len(total_depth)))
			acc_coverage_minDep.append(float(np.count_nonzero(total_depth_minDep)/len(total_depth_minDep)))
	total_coverage = float(np.count_nonzero(total_depth)/len(total_depth))
	total_coverage_minDep = float(np.count_nonzero(total_depth_minDep)/len(total_depth_minDep))


	## plot depth distribution on each sample
	keys = list(mapped_reads.keys())
	for i in range(len(keys)):
		k = keys[i]
		data = depth[k]
		max_d = np.max(data["depth"])
		min_d = np.min(data["depth"])
		median_d = np.median(data["depth"])
		mean_d = np.mean(data["depth"])
		std_d = np.std(data["depth"])
		plt.figure(1, figsize=(20, 8))
		plt.plot(data["pos"], data["depth"])
		plt.axhline(y=minDep, color='r', linestyle='-')
		plt.xlim(0,max_pos)
		plt.ylim(0,max_depth) 
		plt.ylabel("Reads Depth")
		plt.xlabel("Position (in MB)")
		plt.title("Sample " + k +
							"\n Mapped Reads:" + str(mapped_reads[k]) + 
							", Coverage:" + "{:.2%}".format(coverage_list[i]) + 
							", Coverage (depth >= "+ str(minDep)+"):" + "{:.2%}".format(coverage_minDep_list[i]) + 
							"\n Min: " + str(min_d) +
							", Max: " + str(max_d) +
							", Median: " + str(median_d) +
							", Mean: " + str(mean_d) +
							", STD: " + str(std_d))
		plt.savefig(outfile + "/" + k + "_depth.png")
		plt.close()


	## plot depth distribution on all sample
	plt.figure(1, figsize=(20, 8))
	for k in mapped_reads.keys():
		data = depth[k]
		plt.plot(data["pos"], data["depth"])
	plt.ylabel("Reads Depth")
	plt.xlabel("Position (in MB)")
	plt.axhline(y=minDep, color='r', linestyle='-')
	plt.title("Reference Genome " + outfile.split("/")[0] +
						"\n Average Mapped Reads:" + str(np.mean(list(mapped_reads.values()))) + 
						", Average Coverage:" + "{:.2%}".format(np.mean(coverage_list)) + 
						", Average Coverage (depth >= "+ str(minDep)+"):" + "{:.2%}".format(np.mean(coverage_minDep_list)) +
						"\n Total Coverage :" + "{:.2%}".format(total_coverage) +
						", Total Coverage (depth >= "+ str(minDep)+"):" + "{:.2%}".format(total_coverage_minDep))
	plt.xlim(0,max_pos)
	plt.ylim(0,max_depth) 
	plt.savefig(outfile + "/all_depth.png")
	plt.close()


	## plot depth per mapped reads distribution on all sample
	max_depth = 0
	plt.figure(1, figsize=(20, 8))
	for k in mapped_reads.keys():
		data = depth[k]
		plt.plot(data["pos"], data["depth"]/mapped_reads[k])
		if max(data["depth"]/mapped_reads[k]) > max_depth:
				max_depth = max(data["depth"]/mapped_reads[k])
	plt.ylabel("Reads Depth / Mapped Reads")
	plt.xlabel("Position (in MB)")
	plt.title("Reference Genome " + outfile.split("/")[0] +
						"\n Average Mapped Reads:" + str(np.mean(list(mapped_reads.values()))) + 
						", Average Coverage:" + "{:.2%}".format(np.mean(coverage_list))+
						", Total Coverage :" + "{:.2%}".format(total_coverage))
	plt.xlim(0,max_pos)
	plt.ylim(0,max_depth) 
	plt.savefig(outfile + "/all_depth_per_mapped.png")
	plt.close()


	## plot depth (minDep = 100) per mapped reads distribution on all sample 
	max_depth = 0
	plt.figure(1, figsize=(20, 8))
	for k in mapped_reads.keys():
		data = depth[k]
		temp_pos = data["pos"]
		depth_reads = [float(i/mapped_reads[k]) if i >= minDep else 0.0 for i in data["depth"]]
		plt.plot(temp_pos, depth_reads)
		if max(depth_reads) > max_depth:
				max_depth = max(depth_reads)
	plt.ylabel("Reads Depth / Mapped Reads")
	plt.xlabel("Position (in MB)")
	plt.title("Reference Genome " + outfile.split("/")[0] +
						"\n Average Mapped Reads:" + str(np.mean(list(mapped_reads.values()))) + 
						", Average Coverage (depth >= "+ str(minDep)+"):" + "{:.2%}".format(np.mean(coverage_minDep_list))+
						", Total Coverage (depth >= "+ str(minDep)+"):" + "{:.2%}".format(total_coverage_minDep))
	plt.xlim(0,max_pos)
	plt.ylim(0,max_depth) 
	plt.savefig(outfile + "/all_depth_per_mapped_" + str(minDep) + ".png")
	plt.close()


	## plot depth per mapped reads distribution on all sample 
	plt.figure(1, figsize=(20, 8))
	positions = []
	total_mapped_reads = 0
	for k in mapped_reads.keys():
		total_mapped_reads += mapped_reads[k]
		if len(positions) == 0:
			positions = depth[k]["pos"]
	depth_reads = [float(x/total_mapped_reads) for x in total_depth]
	max_d = np.max(depth_reads)
	min_d = np.min(depth_reads)
	median_d = np.median(depth_reads)
	mean_d = np.mean(depth_reads)
	std_d = np.std(depth_reads)
	plt.plot(positions, depth_reads)
	plt.ylabel("Reads Depth / Mapped Reads")
	plt.xlabel("Position (in MB)")
	plt.xlim(0,max_pos)
	plt.ylim(0,max(depth_reads))
	percentage = "{:.2%}".format(np.count_nonzero(total_depth)/len(positions))
	plt.title("Reference Genome " + outfile.split("/")[0] +
						"\n Mapped Reads:" + str(total_mapped_reads) + 
						", Coverage:" + percentage + 
						"\n Min: " + str(min_d) +
						", Max: " + str(max_d) +
						", Median: " + str(median_d) +
						", Mean: " + str(mean_d) +
						", STD: " + str(std_d))
	plt.savefig(outfile + "/total_depth_per_mapped.png")
	plt.close()


	## plot depth (minDep = 100) per mapped reads distribution on all sample 
	plt.figure(1, figsize=(20, 8))
	positions = []
	total_mapped_reads = 0
	for k in mapped_reads.keys():
		total_mapped_reads += mapped_reads[k]
		if len(positions) == 0:
			positions = depth[k]["pos"]
	depth_reads = [float(x/total_mapped_reads) for x in total_depth_minDep]
	max_d = np.max(depth_reads)
	min_d = np.min(depth_reads)
	median_d = np.median(depth_reads)
	mean_d = np.mean(depth_reads)
	std_d = np.std(depth_reads)
	plt.plot(positions, depth_reads)
	plt.ylabel("Reads Depth / Mapped Reads")
	plt.xlabel("Position (in MB)")
	plt.xlim(0,max_pos)
	plt.ylim(0,max(depth_reads))
	percentage = "{:.2%}".format(np.count_nonzero(total_depth_minDep)/len(positions))
	plt.title("Reference Genome " + outfile.split("/")[0] +
						"\n Mapped Reads:" + str(total_mapped_reads) + 
						", Coverage:" + percentage + 
						"\n Min: " + str(min_d) +
						", Max: " + str(max_d) +
						", Median: " + str(median_d) +
						", Mean: " + str(mean_d) +
						", STD: " + str(std_d))
	plt.savefig(outfile + "/total_depth_per_mapped_" + str(minDep) + ".png")
	plt.close()


	## plot accumulated coverage
	plt.figure(1, figsize=(8, 8))
	plt.plot(np.arange(1, len(acc_coverage)+1), acc_coverage, label='min depth = 0')
	plt.plot(np.arange(1, len(acc_coverage_minDep)+1), acc_coverage_minDep, label='min depth = ' + str(minDep))
	plt.legend()
	plt.ylabel("Acc Coverage")
	plt.xlabel("Number of Samples")
	plt.title("Reference Genome " + outfile.split("/")[0] +
			  "\n Accumulated Coverage")
	plt.savefig(outfile + "/acc_coverage.png")
	plt.close()


