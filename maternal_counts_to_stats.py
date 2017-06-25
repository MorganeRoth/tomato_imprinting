#!/usr/local/bin/python
# Python 2.6 
# maternal_counts_to_stats.py
# Morgane Roth, Ana Marcela Florez-Rueda, Mathias Scharmann, PEG Integrative Biology, ETH Zürich 
# January 2016

# example usage:
# python maternal_counts_to_stats.py \
# --rep1_a A1_A13_A25_A10_A44_A34_endosperm_1_counts_per_gene.txt \
# --rep2_a A1_A13_A25_A10_A44_A34_endosperm_2_counts_per_gene.txt \
# --rep3_a A1_A13_A25_A10_A44_A34_endosperm_3_counts_per_gene.txt \
# --rep1_b A1_A13_A25_A10_A44_A34_endosperm_4_counts_per_gene.txt \
# --rep2_b A1_A13_A25_A10_A44_A34_endosperm_5_counts_per_gene.txt \
# --rep3_b A1_A13_A25_A10_A44_A34_endosperm_6_counts_per_gene.txt 
# --output_name

# this script works first with the replicates for each direction of the cross independently
# then it takes the parental signals (obtained/summerized over the 3 reps) of both reciprocals to see if the allele specific expression is due to the parents or just plant-specific
# this is the rationale to find imprinted genes, but the signal is further statistically tested with a chi-square
# this script has 3 important steps:
# - first the median is computed for sites expressed in at least 2 reps for a given direction of the cross
# - then this median is used to find parental bias in each direction of the cross with the rule that the signal should agree in at least 2 libraries/reps
# - by comparing this signal in the reciprocals, we find the "true" imprinted genes.
# the p-value should agree with the results  (the expected proportions are given by the median)

"""
implementing a Chi-squared test 
Do the read counts deviate from the expectation 2M:1P in at least one direction of the cross?
- p-values are FDR-corrected by Benjamini-Hochberg procedure -> q-values
- genes are classified as "not imprinted" if q-value > 0.05

"""
	
import os
import sys
import argparse
import numpy
from scipy.stats import chisquare
import operator

# checks for file and breaks script if not existing
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x


# breaks script if non-UNIX linebreaks in input files
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()
	
# parses command line arguments, breaks script if not plausible
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
		
	parser.add_argument("--rep1_a", required=True, type=extant_file,
		help="name / path of first endopserm file", metavar="FILE")
	parser.add_argument("--rep2_a", required=True, type=extant_file,
		help="name / path of secod endosperm file", metavar="FILE")
	parser.add_argument("--rep3_a", required=True, type=extant_file,
		help="name / path of first endopserm file", metavar="FILE")
	parser.add_argument("--rep1_b", required=True, type=extant_file,
		help="name / path of secod endosperm file", metavar="FILE")
	parser.add_argument("--rep2_b", required=True, type=extant_file,
		help="name / path of first endopserm file", metavar="FILE")
	parser.add_argument("--rep3_b", required=True, type=extant_file,
		help="name / path of secod endosperm file", metavar="FILE")
	parser.add_argument("--output_name", required=True)
	
	args = parser.parse_args()
	
	linebreak_check(args.rep1_a)
	linebreak_check(args.rep2_a)
	linebreak_check(args.rep3_a)
	linebreak_check(args.rep1_b)
	linebreak_check(args.rep2_b)
	linebreak_check(args.rep3_b)

	# finish
	return args
	
#############

# this function will test if the gene is expressed in at least 2 reps and returns a ratio MAT/TOT counts based by summing all counts
# this will not be used for individual genes, but to get the median and so the thresholds for imprinted genes
 
def add_counts (a_data, b_data, c_data):
	
	all_genes = set(c_data.keys()).union(set(a_data.keys()).union(set(b_data.keys())))
# we add all counts with a union	
	print "rep1 genes found:	", len(set(a_data.keys()))
	print "rep2 genes found:	", len(set(b_data.keys()))
	print "rep3 genes found:	", len(set(c_data.keys()))
	print "genes found in intersection of 3 reps:	", len(set(a_data.keys()).intersection(set(b_data.keys())).intersection(set(c_data.keys())))
	print "total genes found:	", len(all_genes)
# we check that the count is done at least twice	
	sumcounts = {}
	for gene in all_genes:
		try:
			a_mat = float(a_data[gene][0])
			a_tot = float(a_data[gene][1])
			rep_a=1			
		except KeyError:
			a_mat = 0.0
			a_tot = 0.0
			rep_a=0
		try:
			b_mat = float(b_data[gene][0])
			b_tot = float(b_data[gene][1])
			rep_b=1
		except KeyError:
			b_mat = 0.0
			b_tot = 0.0
			rep_b=0
		try:
			c_mat = float(c_data[gene][0])
			c_tot = float(c_data[gene][1])
			rep_c=1
		except KeyError:
			c_mat = 0.0
			c_tot = 0.0
			rep_c=0
# when a gene is expressed at least in 2 libraries, we count the maternal reads, the total reads, and the number of libraries where it is found
		if not rep_a + rep_b + rep_c == 1:
			ratio= (a_mat + b_mat + c_mat)/ (a_tot + b_tot + c_tot)
			sumcounts[gene] = [a_mat + b_mat + c_mat, a_tot + b_tot + c_tot, ratio, rep_a + rep_b + rep_c]
		rep_a=0
		rep_b=0
		rep_c=0
		
	return sumcounts	
	
# in the sumcounts element, the third element (ratio), will be used for further calculations, but the other elements will be output in the "added" file for drawing graphs

def write_output (outdata, header, outfile):

	outlines = ["\t".join(header)]
	for gene in sorted(outdata.keys()):
		outlines.append("\t".join([gene, str(outdata[gene][0]), str(outdata[gene][1]),str(outdata[gene][2]), str(outdata[gene][3])]))
	with open(outfile, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))
		
# this function prints all ratios in a list and calculates the median
def read_ratios (added_data):
	all_genes = set(added_data.keys())
	mylist=[]
	count=0
	for gene in all_genes:
		ratio=added_data[gene][2]
		mylist.append(ratio)
		count +=1 
	median=numpy.median(mylist)
	print "the median of the maternal proportion is", median
	print "observed in", count, "genes and non-annotated sites"
	return median

# one median is calculated for each cross (3reps), and we compute the thresholds used for 
# parent specific expression detection (not called imprinting yet) in one direction of the cross
	
def compute_threshold (median):
	pat_prop=1-float(median)
	pat_increase_by_two=pat_prop*2 
	pat_decrease_by_two=pat_prop/2 
	threshold_meg=1-pat_decrease_by_two
	threshold_peg=float(median)/2
	threshold_MEG=1-(pat_decrease_by_two/2)
	threshold_PEG=threshold_peg/2
	print "pat_prop" , pat_prop
	print "threshold_meg", threshold_meg, "threshold_peg", threshold_peg, "threshold_MEG", threshold_MEG, "threshold_PEG", threshold_PEG
	
	return threshold_meg, threshold_peg, threshold_MEG, threshold_PEG

# we import the read counts for all reps

def read_to_dict (file):	
	data = {}
	with open(file, "r") as INFILE:
		header = INFILE.readline().strip("\n").split("\t")
# 		print header
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			data[str(fields[0])] = fields[1:]
	return data, header

# we compare the ratio for each gene/location to the computed threshold and assign the status
# NI = non imprinted
# meg, peg = maternally vs. paternally expressed
# MEG, PEG = strongly biased expression towards one parent
# we have to run this function twice (for each direction of the cross)

def compare_ratio (threshold_meg, threshold_peg, threshold_MEG, threshold_PEG, rep1_data, rep2_data, rep3_data):
	all_genes = set(rep1_data.keys()).union(set(rep2_data.keys()).union(set(rep3_data.keys())))
	summary={}
	for gene in all_genes:
		tot_meg=0
		tot_peg=0
		tot_MEG=0
		tot_PEG=0		
		rep1_mat=0
		rep1_tot=0
		rep2_mat=0
		rep2_tot=0
		rep3_mat=0
		rep3_tot=0
		count_1_meg=0
		count_2_meg=0
		count_3_meg=0
		count_1_peg=0
		count_2_peg=0
		count_3_peg=0
		count_1_MEG=0
		count_2_MEG=0
		count_3_MEG=0
		count_1_PEG=0
		count_2_PEG=0
		count_3_PEG=0
		status=0
		try:
			rep1_mat = float(rep1_data[gene][0])
			rep1_tot = float(rep1_data[gene][1])
			ratio1=float(rep1_mat/rep1_tot)
			if (float(ratio1) >= float(threshold_meg)):
				count_1_meg=1
				if (float(ratio1) >= float(threshold_MEG)):
					count_1_MEG=1
			elif (float(ratio1) <= float(threshold_peg)):
				count_1_peg=1
				if (float(ratio1)<= float(threshold_PEG)):
					count_1_PEG=1
			elif	( float(threshold_peg) <= float(ratio1) <= float(threshold_meg) ):
					count_1_meg=0
					count_1_peg=0
					count_1_MEG=0
					count_1_PEG=0			
		except KeyError:
			rep1_tot == 0
			ratio1="NA"
		try:	
			rep2_mat = float(rep2_data[gene][0])
			rep2_tot = float(rep2_data[gene][1])
			ratio2=float(rep2_mat/rep2_tot)
			if ratio2 >= threshold_meg:
				count_2_meg=1
				if ratio2 >=threshold_MEG:
					count_2_MEG=1
			elif ratio2 <= threshold_peg:
				count_2_peg=1
				if ratio2 <=threshold_PEG:
					count_2_PEG=1	
			elif  ( float(threshold_peg) <= float(ratio2) <= float(threshold_meg) ):
					count_2_meg=0
					count_2_peg=0
					count_2_MEG=0
					count_2_PEG=0							
		except KeyError:
			rep2_tot == 0
			ratio2="NA"
		try:
			rep3_mat = float(rep3_data[gene][0])
			rep3_tot = float(rep3_data[gene][1])
			ratio3=float(rep3_mat/rep3_tot)
			if ratio3 >= threshold_meg:
				count_3_meg=1
				if ratio3 >=threshold_MEG:
					count_3_MEG=1
			elif ratio3 <= threshold_peg:
				count_3_peg=1
				if ratio3 <=threshold_PEG:
					count_3_PEG=1
			elif  ( float(threshold_peg) <=float(ratio3) <= float(threshold_meg) ):
					count_3_meg=0
					count_3_peg=0
					count_3_MEG=0
					count_3_PEG=0	
		except KeyError:
			rep3_tot == 0
			ratio3="NA"
		tot_meg=count_1_meg + count_2_meg + count_3_meg
		tot_peg=count_1_peg + count_2_peg + count_3_peg
		tot_MEG=count_1_MEG + count_2_MEG + count_3_MEG
		tot_PEG=count_1_PEG + count_2_PEG + count_3_PEG
		if (tot_meg in range(2)) and (tot_peg in range(2)) :
			status= "NI"
		elif tot_meg >= 2 :
			status = "meg"
			if tot_MEG >= 2 :
				status = "MEG"
		elif tot_peg >= 2:
			status = "peg"
			if tot_PEG >= 2 :
				status= "PEG"
		summary[gene] = [ratio1, ratio2, ratio3, tot_meg, tot_peg, tot_MEG, tot_PEG, status]
		tot_meg=0
		tot_peg=0
		tot_MEG=0
		tot_PEG=0		
		rep1_mat=0
		rep1_tot=0
		rep2_mat=0
		rep2_tot=0
		rep3_mat=0
		rep3_tot=0
		count_1_meg=0
		count_2_meg=0
		count_3_meg=0
		count_1_peg=0
		count_2_peg=0
		count_3_peg=0
		count_1_MEG=0
		count_2_MEG=0
		count_3_MEG=0
		count_1_PEG=0
		count_2_PEG=0
		count_3_PEG=0
		
	return summary

# now we make a file with data for both directions of the cross

def combine_signals_recip (summary_a, summary_b):
	intersection_genes = set(summary_a.keys()).intersection(set(summary_b.keys()))
	summary_ab ={}
	for gene in intersection_genes:

		for i in range(8):
			summary_ab[gene]=summary_a[gene]
		for i in range(8):
			summary_ab[gene].append(summary_b[gene][i])
	return summary_ab
		
# in the main, we combine all columns from summary_direction1 and summary_direction2
# we use this file to compare the status of each gene in both directions and conclude for imprinted genes (columns 7 and 15)

def compare_status (combined_reps_data):
	# take only the common genes
	intersection_genes = set(combined_reps_data.keys())
	overlap={}
	for gene in intersection_genes:
		status_ab= "null"
		status_a = str(combined_reps_data[gene][7])
		status_b = str(combined_reps_data[gene][15])
		if status_a == status_b:
			status_ab = status_a
		elif (status_a == "meg" and status_b == "MEG") or (status_a == "MEG" and status_b == "meg"):
			status_ab = "meg"
		elif  (status_a == "peg" and status_b == "PEG") or (status_a == "PEG" and status_b == "peg"):
			status_ab = "peg"
		elif  (status_a == "peg" and status_b == "MEG") or (status_a == "peg" and status_b == "meg"):
			status_ab = "plant_specific"
		elif  (status_a == "NI") or (status_b == "NI"):
			status_ab = "NI"
		else :
			status_ab = "plant_specific"	
		overlap[gene] = status_ab
		
	return overlap

# make a chisquare test : expected values are based on the median

def do_chisq (added_data_a, added_data_b, median_a, median_b ):
	
	all_genes = set(added_data_a.keys()).union(set(added_data_b.keys()))
	
	print "endosperm a genes found:	", len(set(added_data_a.keys()))
	print "endosperm b (reciprocal of a) genes found:	", len(set(added_data_b.keys()))
	print "overlap of expressed genes in the endosperms of the reciprocal cross:	", len(set(added_data_a.keys()).intersection(set(added_data_b.keys())))
		
	outdict = {}
	for gene in all_genes:
		try:
			a_ratio = added_data_a[gene][2]
		except KeyError:
			continue
			# this ensures that genes which have not been 
			# seen in BOTH reciprocal crosses' endosperms are NOT put into the output
		try:
			b_ratio = added_data_b[gene][2]
		except KeyError:
			continue 
			# this ensures that genes which have not been 
			# seen in BOTH reciprocal crosses' endosperms are NOT put into the output
		
		# Chi sqaured test on the read counts vs direction of the cross
		x= added_data_a[gene][1] * median_a
		ma= round(added_data_a[gene][1] * median_a)
		pa= round(added_data_a[gene][1] * (1-median_a))
		mb= round(added_data_b[gene][1] * median_b)
		pb= round(added_data_b[gene][1] * (1-median_b))		
		obs_counts = [added_data_a[gene][0],added_data_a[gene][1]-added_data_a[gene][0],added_data_b[gene][0],added_data_b[gene][1]-added_data_b[gene][0]]
		exp_counts = [ma, pa, mb, pb]
		test_result = chisquare(f_obs = obs_counts, f_exp = numpy.array(exp_counts)) 
		p_value = test_result[1]
		outdict[gene] = [a_ratio, b_ratio, p_value]
	return outdict	

# Now that we have the p-values, we compute the Benjamini-Hochberg correction >> q-values

def bh_qvalues(p_values):
	"""
	Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

	This function implements an algorithm equivalent to L{bh_rejected} but
	yields a list of 'adjusted p-values', allowing for rejection decisions
	based on any given threshold.

	@type pv: list
	@param pv: p-values from a multiple statistical test

	@rtype: list
	@return: adjusted p-values to be compared directly with the desired FDR
	  level

	# Copyright (C) 2008 Simone Leo - CRS4. All Rights Reserved.
	# 
	# Permission to use, copy, modify, and distribute this software and its
	# documentation for educational, research, and not-for-profit purposes, without
	# fee and without a signed licensing agreement, is hereby granted, provided
	# that the above copyright notice, this paragraph and the following two
	# paragraphs appear in all copies, modifications, and distributions. Contact
	# CRS4, Parco Scientifico e Tecnologico, Edificio 1, 09010 PULA (CA - Italy),
	# +3907092501 for commercial licensing opportunities.
	# 
	# IN NO EVENT SHALL CRS4 BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
	# INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
	# THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF CRS4 HAS BEEN ADVISED
	# OF THE POSSIBILITY OF SUCH DAMAGE.
	# 
	# CRS4 SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
	# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
	# PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
	# HEREUNDER IS PROVIDED "AS IS". CRS4 HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
	# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
	
	this function was downloaded from:
	http://pydoc.net/Python/ngslib/1.1.4/ngslib.fdr/
	05 Sep 2015

	"""
	if not p_values:
		return []
	m = len(p_values)
	args, p_values = zip(*sorted(enumerate(p_values), None, operator.itemgetter(1)))
	if p_values[0] < 0 or p_values[-1] > 1:
		raise ValueError("p-values must be between 0 and 1")
	q_values = m * [0]
	mincoeff = p_values[-1]
	q_values[args[-1]] = mincoeff
	for j in xrange(m-2, -1, -1):
		coeff = m*p_values[j]/float(j+1)
		if coeff < mincoeff:
			mincoeff = coeff
		q_values[args[j]] = mincoeff
	return q_values

def do_FDR (p_value_dict, summary_ab, overlap ):
	
	genes = p_value_dict.keys()
	
	p_vals = [p_value_dict[x][2] for x in genes]

# p_values are transformed in q-value by a Benjamini Hochberg correction (bh)
	qvals = bh_qvalues(p_vals)
	
	outdict = {}
	for i, gene in enumerate(genes):
		outdict[gene] = summary_ab[gene]
		outdict[gene].append(overlap[gene])
		outdict[gene].append(p_value_dict[gene][0])
		outdict[gene].append(p_value_dict[gene][1])
		outdict[gene].append(p_value_dict[gene][2])
		outdict[gene].append(qvals[i])
		
	return outdict
	

def write_output2 (outdata, header, outfile):

	outlines = ["\t".join(header)]
	for gene in sorted(outdata.keys()):
		list=[]
		list=[gene]
		for i in range (21):
			list.append(str(outdata[gene][i]))
		outlines.append("\t".join(list))
	with open(outfile, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))
		
		
######### MAIN ########

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
		
	args = 	get_commandline_arguments ()

	
	rep1_a_data, rep1_a_header = read_to_dict (args.rep1_a)
	rep2_a_data, rep2_a_header = read_to_dict (args.rep2_a)
	rep3_a_data, rep3_a_header = read_to_dict (args.rep3_a)
	rep1_b_data, rep1_b_header = read_to_dict (args.rep1_b)
	rep2_b_data, rep2_b_header = read_to_dict (args.rep2_b)
	rep3_b_data, rep3_b_header = read_to_dict (args.rep3_b)
	
	added_data_a = add_counts (rep1_a_data, rep2_a_data, rep3_a_data)
	added_data_b = add_counts (rep1_b_data, rep2_b_data, rep3_b_data)
	

	write_output (added_data_a, rep1_a_header, "added_counts.{0}.{1}.{2}.txt".format(args.rep1_a, args.rep2_a, args.rep3_a))
	write_output (added_data_b, rep1_b_header, "added_counts.{0}.{1}.{2}.txt".format(args.rep1_b, args.rep2_b, args.rep3_b))


	median_a = read_ratios (added_data_a)
	median_b = read_ratios (added_data_b)
	
	threshold_meg_a, threshold_peg_a, threshold_MEG_a, threshold_PEG_a = compute_threshold (median_a)
	threshold_meg_b, threshold_peg_b, threshold_MEG_b, threshold_PEG_b = compute_threshold (median_b)

	summary_a = compare_ratio (threshold_meg_a, threshold_peg_a, threshold_MEG_a, threshold_PEG_a, rep1_a_data, rep2_a_data, rep3_a_data)
	summary_b = compare_ratio (threshold_meg_b, threshold_peg_b, threshold_MEG_b, threshold_PEG_b, rep1_b_data, rep2_b_data, rep3_b_data)
	summary_ab = combine_signals_recip (summary_a, summary_b)
	
	overlap = compare_status (summary_ab)
	
	p_value_dict = do_chisq (added_data_a, added_data_b, median_a, median_b )
		
	q_value_dict = do_FDR (p_value_dict, summary_ab, overlap)
	header = ["gene", "rep1_a", "rep2_a", "rep3_a", "meg_a", "peg_a", "MEG_a", "PEG_a", "status_a", "rep1_b", "rep2_b", "rep3_b", "meg_b", "peg_b", "MEG_b", "PEG_b", "status_b", "status_reciprocal", "a_added_counts_ratio", "b_added_counts_ratio", "chisquare_p_value", "q_value"]
	
	write_output2 (q_value_dict, header, args.output_name)

	print "Done!"
	
	
if __name__ == "__main__":
	main(sys.argv[1:])
