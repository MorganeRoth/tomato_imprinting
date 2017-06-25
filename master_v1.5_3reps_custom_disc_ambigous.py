#!/usr/local/bin/python
# Python 2.6
# master_v1.5_3reps_custom_disc_ambigous.py
# Morgane Roth, Ana Marcela Florez-Rueda, Mathias Scharmann, PEG Integrative Biology, ETH Zürich 
# 01 October 2015

# example usage:
# python /gdc_home4/morgane/softwares/master_v1.5_3reps_custom_disc_ambigous.py \
# --parent_rc LA2185A_LA4329B_minA2_mincov5_rc \
# --endo_rc A6_A18_A30_A4_A16_A28_minA1_mincov2_rc \
# --parent_annot LA2185A_LA4329B_minA2_mincov5_rc.annotated.txt.cleaned.txt  \
# --o A6_A18_A30_A4_A16_A28 \
# --both_hom_tresh 0.05 \
# --het_het_tresh 0.3 \
# --het_hom_tresh 0.05
	
import os
import sys
import argparse
import numpy
 
class Printer():
	"""
	Print things to stdout on one line dynamically
	"""
	def __init__(self,data):
 
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

# check files exist in the folder
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

# check non-UNIX linebreaks 
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()
	
# check command line arguments 
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
		
	parser.add_argument("--parent_rc", required=True, type=extant_file,
		help="name / path of parent rc file", metavar="FILE")
	parser.add_argument("--endo_rc", required=True, type=extant_file,
		help="name / path of endosperm rc file", metavar="FILE")
	parser.add_argument("--parent_annot", required=True, type=extant_file,
		help="name / path of parent annotation file, prouced by perl script then cleaned by python script", metavar="FILE")
	parser.add_argument("--outfile", "--o", required=True, 
		help="name of outfile")
	parser.add_argument("--both_hom_tresh", required=True, 
		help="""both_hom_tresh = for sites that are fixed for different alleles in both parents ("homozygous sites"): maximum % of reads of the minor allele for a site to be called homozygous""", type=float)
	parser.add_argument("--het_het_tresh", required=True, 
		help="""het_het_tresh = for sites where one parent is het and the other hom ("heterozygous sites"): minimum % of reads of the minor allele for a site to be called heterozygous in the heterozygote parent""", type=float)
	parser.add_argument("--het_hom_tresh", required=True, 
		help="""het_hom_tresh = for sites where one parent is het and the other hom ("heterozygous sites"): maximum % of reads of the minor allele for a site to be called homozygous in the homozygote parent""", type=float)
	
	args = parser.parse_args()
	linebreak_check(args.parent_rc)
	linebreak_check(args.endo_rc)
	linebreak_check(args.parent_annot)
	
	# finish
	return args
	
#############

## read data in usable format
def read_to_2dimlist (file):
	
	data = []
	count = 0
	with open(file, "r") as INFILE:
		header = INFILE.readline().strip("\n").split("\t")
		for line in INFILE:
			data.append(line.strip("\n").split("\t"))
			count += 1
	print "infile sites:	", count
		
	return data, header
	
def read_to_dict (file):
	
	data = {}
	with open(file, "r") as INFILE:
		header = INFILE.readline().strip("\n").split("\t")
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			linedict = {}
			for field, head in zip(fields[2:], header[2:]):
				linedict[head] = field
			data["\t".join([fields[0], fields[1]])] = linedict
	
	return data, header

###########

# define different categories of sites and identify them in the parental comparisons

def filter_sites (rc_data, both_hom_tresh, het_het_tresh, het_hom_tresh):
	"""
	both_hom_tresh = for sites that are fixed for different alleles in both parents ("homozygous sites"): maximum % of reads of the minor allele for a site to be called homozygous
	het_het_tresh = for sites where one parent is het and the other hom ("heterozygous sites"): minimum % of reads of the minor allele for a site to be called heterozygous in the heterozygote parent
	het_hom_tresh = for sites where one parent is het and the other hom ("heterozygous sites"): maximum % of reads of the minor allele for a site to be called homozygous in the homozygote parent
	
	"""
	# check that ratio of minor to major allele is less than 5%
	both_hom = []
	plant_1_heto = []
	plant_2_heto = []
	homcount = 0
	hetcount = 0
	for line in rc_data:
		#print line
		# discard lines with "N" in major alleles
		if "N" not in line[7]:
			rc1 = [float(x) for x in line[9].split("/")]
			rc2 = [float(x) for x in line[10].split("/")]
			min_maj_ratio1 = (rc1[1]-rc1[0])/rc1[0]
			min_maj_ratio2 = (rc2[1]-rc2[0])/rc2[0]
			
			# discard lines where major alleles are identical	
			if line[7][0] != line[7][1]:
				# check that ratio of minor to major allele is less than X% (both_hom_treshold) in both libraries
				if min_maj_ratio1 < both_hom_tresh and min_maj_ratio2 < both_hom_tresh:
					both_hom.append(line)
					homcount += 1
			# check that ratio of minor to major allele is < het_hom_tresh in one plant but > het_het_tresh in the other
			if min_maj_ratio1 > het_het_tresh and min_maj_ratio2 < het_hom_tresh:
				plant_1_heto.append(line)
				hetcount += 1
			elif min_maj_ratio1 < het_hom_tresh and min_maj_ratio2 > het_het_tresh:
				plant_2_heto.append(line)
				hetcount += 1
	
	print "sites where both plants are homozygous for different alleles:	", homcount
	print "sites where 1 of 2 plants is heterozygous:	", hetcount 
	return both_hom, plant_1_heto, plant_2_heto

# return sites	
def get_the_sites (inlist):
	
	sites_set = set()
	for line in inlist:
		siteID = "\t".join([line[0], line[1]])
		sites_set.add(siteID)
	
	return sites_set

# count sites found sites in the 6 endosperm files
def filter_endo (target_sites, endodata):

	count = 0
	outendo = []
	for line in endodata:
		siteID = "\t".join([line[0], line[1]])
		if siteID in target_sites:
			outendo.append(line)
			count += 1
	print "parents sites found in endo:	", count
	return outendo
# identify sites to make counts
def find_mat_pat_alleles (sites, parent_data_dict, mothercolumn, fathercolumn):	
	endodict = {}
	for siteID in sorted(sites):
		endodict[siteID] = {"mat":parent_data_dict[siteID][7][mothercolumn], "pat":parent_data_dict[siteID][7][fathercolumn]}
	
	return endodict

# identify the "fixed" and "discriminant" bases for hom-het sites (for example AA-AC; A is fixed, C is discriminant)	
def find_discr_fixed_bases (plant_1_heto_dict, plant_2_heto_dict):	
	# fixed_base_position_in_col7 is 0 for plant_2_heto and 1 for plant_1_heto 
	discr_fixed_dict = {}
	for siteID in sorted(plant_1_heto_dict.keys()):
		# print plant_1_heto_dict[siteID]
		fixed_base = plant_1_heto_dict[siteID][7][1]
		other_base_pos = 0	
		# find discr_base in either col7 or col8 of parent_data_dict
		# if major alleles (column 7) are the same, the discriminant allele is to be found in the minor alleles, column 8
		# the coverage is to be found in either column 9, or in column 11 for plant1 (col 10 and 12 for plant2)
		if plant_1_heto_dict[siteID][7][1] == plant_1_heto_dict[siteID][7][other_base_pos]:
			discr_base = plant_1_heto_dict[siteID][8][other_base_pos]
			discr_cov= plant_1_heto_dict[siteID][11].split("/")[0]
		else:
			discr_base = plant_1_heto_dict[siteID][7][other_base_pos]
			discr_cov= plant_1_heto_dict[siteID][9].split("/")[0]
		discr_fixed_dict[siteID] = {"discr":discr_base, "fixed":fixed_base, "discr_cov":discr_cov}
	
	for siteID in sorted(plant_2_heto_dict.keys()):
		fixed_base = plant_2_heto_dict[siteID][7][0]
		other_base_pos = 1	
		# find discr_base in either col7 or col8 of parent_data_dict
		if plant_2_heto_dict[siteID][7][0] == plant_2_heto_dict[siteID][7][other_base_pos]:
			discr_base = plant_2_heto_dict[siteID][8][other_base_pos]
			discr_cov= plant_2_heto_dict[siteID][12].split("/")[0]
		else:
			discr_base = plant_2_heto_dict[siteID][7][other_base_pos]
			discr_cov= plant_2_heto_dict[siteID][10].split("/")[0]
		
		discr_fixed_dict[siteID] = {"discr":discr_base, "fixed":fixed_base, "discr_cov":discr_cov}
	
	return discr_fixed_dict		 

def make_dict_from_2dimlist(inlist):
	
	outdict = {}
	for line in inlist:
		siteID = "\t".join([line[0], line[1]])
		outdict[siteID] = line
	
	return outdict

def get_endo_counts_mat_pat (endosperm, target_sites, endo_data_dict, endo_mat_pat_dict, endo_info):
	
	outcounts = {}
	for siteID in endo_data_dict.keys():
		if siteID in target_sites:
			maj_all = endo_data_dict[siteID][7][int(endo_info[endosperm]["endocolumn"])-9]
			min_all = endo_data_dict[siteID][8][int(endo_info[endosperm]["endocolumn"])-9]
			if maj_all == endo_mat_pat_dict[siteID]["mat"]:
				outcounts[siteID] = endo_data_dict[siteID][int(endo_info[endosperm]["endocolumn"])]
				#print "MAJOR maternal"
			elif min_all == endo_mat_pat_dict[siteID]["mat"]:
				outcounts[siteID] = endo_data_dict[siteID][int(endo_info[endosperm]["endocolumn"]+6)]
				#print "minor maternal"
				
	return outcounts
	
def get_mother_zyg (endosperm, endo_info, p1_het_sites, p2_het_sites):
	
	if endo_info[endosperm]["mothercolumn"] == 9:
		motherhetset = p1_het_sites
		motherhomset = p2_het_sites
	elif endo_info[endosperm]["mothercolumn"] == 10:
		motherhetset = p2_het_sites
		motherhomset = p1_het_sites
	
	return motherhetset, motherhomset	

# get counts from hom-het sites; if a discriminant site is not found in the endosperm,
# it could have been a sequencing error in the parental bud transcriptome 
# >>>> thus for these cases, filter sites with <10 reads coverage in parents for the discriminant allele

def get_n_transform_het_counts (endo_data_dict, motherhetset, motherhomset, endo_info, endosperm, discr_fixed_dict):
	
	all_het_sites = motherhetset.union(motherhomset)
	
	outdict = {}
	for siteID in endo_data_dict.keys():
		
		if siteID in all_het_sites: # this filter before all the rest to avoid processing of totally unused sites
			maj_all = endo_data_dict[siteID][7][int(endo_info[endosperm]["endocolumn"])-9]
			min_all = endo_data_dict[siteID][8][int(endo_info[endosperm]["endocolumn"])-9]
			if maj_all == discr_fixed_dict[siteID]["discr"]:
				discr_to_total = endo_data_dict[siteID][int(endo_info[endosperm]["endocolumn"])]
			elif min_all == discr_fixed_dict[siteID]["discr"]:
				discr_to_total = endo_data_dict[siteID][int(endo_info[endosperm]["endocolumn"]+6)]
			elif maj_all == "N": 
			# this is a filtering step to remove sites with zero reads: major allele can never be N
			# if there are any reads at all for this site; minor allele CAN be N!
#				print "jumped to next site because maj alleles was N / zero reads for site"
				continue	
			elif min_all == "N" and float(discr_fixed_dict[siteID]["discr_cov"])<10:  
			# this step removed sites where we are not sure about the discriminant in the parents,
			# because it is not found in the endosperm
			# if the allele is not observed in the endosperm, it has to have enough coverage/confidence in the parents
				continue
			elif min_all == "N" and float(discr_fixed_dict[siteID]["discr_cov"])>=10: 
# 			and discrimant_coverage >= 10:
			# if the minor allele is N, it has zero reads. We only get to this if-level if the major 
			# allele is not the discriminant base AND the minor is N because it has no reads in the endosperm.
            # This in turn implies that the discriminant base has zero reads!
            # Then we have to check that the coverage of the discriminant in the parents has enough coverage
            # because if it is not found in the endosperm, it could be  sequencing error (threshold of 10 counts)
# 				print discr_fixed_dict[siteID]
				total = endo_data_dict[siteID][int(endo_info[endosperm]["endocolumn"])].split("/")[1]
				discr_to_total = "0/{0}".format(total)
#				print "set discr count to zero:	", discr_to_total
				
			if siteID in motherhetset:
			# if mother is heterozygous for this site, then the discriminant base is from the mother
			# -> get count of discriminant base and transform proportion as prop*2	
#				print siteID, "mother het"
				newcount = int(discr_to_total.split("/")[0])*2
				t = int(discr_to_total.split("/")[1])
				if newcount > t: # cap the transformed count at the total count in case it exceeds the total count (to be reasonable)!
					outdict[siteID] = str(t) + "/" + str(t)
				else:
					outdict[siteID] = str(newcount) + "/" + str(t)
#				print "oldcount:	", int(discr_to_total.split("/")[0]), "newcount:	", newcount
						
			elif siteID in motherhomset:
			# if mother is homozygous for this site, then the discriminant base is from the father 
			# -> get count of discriminant base and transform proportion as 1-(prop*2)
#				print siteID, "mother HOM"
				t = discr_to_total.split("/")[1]
				prop = (float(discr_to_total.split("/")[0]))/(float(t))
				newprop = 1-(prop*2)
				if newprop < 0: # if negative, cap at 0 to be reasonable!
					outdict[siteID] = str(0.0) + "/" + t
					fufu = outdict[siteID]
#					print "negative found"
				else:
					outdict[siteID] = str(newprop*float(t)) + "/" + t
#				print "oldcount:	", discr_to_total.split("/")[0], "newcount:	", str(newprop*float(t)), "	tis one was 1-(prop*2)"
				
				
	return outdict	

# write output

def write_output (outdata, rc_header, outfile):
	
	outlines = ["\t".join(rc_header)]
	for line in outdata:
		outlines.append("\t".join(line))
	
	with open(outfile, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def write_dict_to_file (outdict, outfilename):
	
	outheader = ["##chr", "pos", "maternal/total"]
	outlines = ["\t".join(outheader)]
	for siteID in sorted(outdict.keys()):
		chr = siteID.split("\t")[0]
		pos = siteID.split("\t")[1]
		mt = outdict[siteID]
		outlines.append("\t".join([chr, pos, mt]))
	
	with open (outfilename, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def read_gtf(gtf_file):

	gtf = {}
	with open(gtf_file, "r") as INFILE:
		gene_ID = 0
		for line in INFILE:
			gene_ID += 1
			fields = line.strip("\n").split("\t")
			chr = fields[0]
			start = int(fields[3])
			end = int(fields[4])
			if chr not in set(gtf.keys()):
				gtf[chr] = {gene_ID:{"data":fields, "start":start, "end":end}}
			else:
				gtf[chr][gene_ID] = {"data":fields, "start":start, "end":end}
	
	print "gtf file parsed"
	print len(gtf)
	
	return gtf

# use annotation to gather counts per gene

def get_annotation (sites, parent_annot):

	out = {} 
	counter = 0
	totlen = len(sites)
	print "annotating sites to genes, percent progress: "
	for site in sites:
		chr = site.split("\t")[0]
		pos = int(site.split("\t")[1])
		try:
			gene = parent_annot[site]
			out[site] = gene
		except KeyError:
			out[site] = "not annotated in the parents"
			pass
		
		counter += 1
		Printer(float(counter)/float(totlen)*100)
		
	print "got annotated endopserm from parents annotation"
	return out

def count_per_gene (all_sites, endo_sites_annotated):
	
	outcounts = {}
	for gene in endo_sites_annotated.values():
		outcounts[gene] = [0,0,0,""]
	
	for siteID, count in all_sites.items():
		gene = endo_sites_annotated[siteID]
		mat_count = float(count.split("/")[0])
		tot_count = float(count.split("/")[1])
		outcounts[gene][0] += mat_count 
		outcounts[gene][1] += tot_count
		outcounts[gene][2] += 1
		outcounts[gene][3] = outcounts[gene][3] + "<>" + siteID
		
	print "counts per gene finished"
	return outcounts

def write_genecounts_to_file (outdict, outfilename):
	
	outheader = ["gene_name", "maternal", "total", "number_of_sites", "sites"]
	outlines = ["\t".join(outheader)]
	for gene in sorted(outdict.keys()):
		mat = str(outdict[gene][0])
		tot = str(outdict[gene][1])
		nsites = str(outdict[gene][2])
		sites = str(outdict[gene][3])
		outlines.append("\t".join([gene, mat, tot, nsites, sites]))
	
	with open (outfilename, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def read_aprent_annot_to_dict (file):
	
	data = {}
	with open(file, "r") as INFILE:
		INFILE.readline # discard header
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			data["\t".join([fields[0], fields[1]])] = fields[2]	
	return data

	
########################

######### MAIN ########

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
		
	args = 	get_commandline_arguments ()
	
	print args
	
	rc_data, rc_header = read_to_2dimlist (args.parent_rc)
	
	both_hom, plant_1_heto, plant_2_heto = filter_sites (rc_data, args.both_hom_tresh, args.het_het_tresh, args.het_hom_tresh)
	
	endo_data, endo_header = read_to_2dimlist (args.endo_rc)
	
	both_hom_sites = get_the_sites (both_hom)
	p1_het_sites = get_the_sites (plant_1_heto)
	p2_het_sites = get_the_sites (plant_2_heto)
	
	homendo =  filter_endo(both_hom_sites, endo_data)
	
	#######
	
	endo_info = {"1":{"endocolumn":9,"mothercolumn":9, "fathercolumn":10, "col7mother":0, "col7father":1},
	"2":{"endocolumn":10,"mothercolumn":9, "fathercolumn":10, "col7mother":0, "col7father":1},
	"3":{"endocolumn":11,"mothercolumn":9, "fathercolumn":10, "col7mother":0, "col7father":1},
	"4":{"endocolumn":12,"mothercolumn":10, "fathercolumn":9, "col7mother":1, "col7father":0},
	"5":{"endocolumn":13,"mothercolumn":10, "fathercolumn":9, "col7mother":1, "col7father":0},
	"6":{"endocolumn":14,"mothercolumn":10, "fathercolumn":9, "col7mother":1, "col7father":0},
	}
	
	parent_data_dict = make_dict_from_2dimlist(rc_data)
	plant_1_heto_dict = make_dict_from_2dimlist(plant_1_heto)
	plant_2_heto_dict = make_dict_from_2dimlist(plant_2_heto)
	endo_data_dict = make_dict_from_2dimlist(endo_data)
	
	discr_fixed_dict = find_discr_fixed_bases (plant_1_heto_dict, plant_2_heto_dict)
	
	parent_sites_annotated = read_aprent_annot_to_dict (args.parent_annot)
	
	for endosperm in [1,2,3,4,5,6]:
	
		mothercolumn = endo_info[str(endosperm)]["col7mother"] # the position of the mother in column 7 of the parent_rc file
		fathercolumn = endo_info[str(endosperm)]["col7father"] # the position of the father in column 7 of the parent_rc file
		
		# first the sites where both parents are homozygous:
		endo_mat_pat_dict = find_mat_pat_alleles (both_hom_sites, parent_data_dict, mothercolumn, fathercolumn)
		endocounts_homsites = get_endo_counts_mat_pat (str(endosperm), both_hom_sites, endo_data_dict, endo_mat_pat_dict, endo_info)
		
		# now the sites where one parent is heterozygous (has discriminant base) and the other is homozygous and the coverage of the discriminant:
		motherhetset, motherhomset = get_mother_zyg (str(endosperm), endo_info, p1_het_sites, p2_het_sites)
		
		hetsites_transf_counts = get_n_transform_het_counts (endo_data_dict, motherhetset, motherhomset, endo_info, str(endosperm), discr_fixed_dict)
		
		# finish:
		
		all_sites = dict(endocounts_homsites.items() + hetsites_transf_counts.items())
		outfilename = "{0}.endosperm_{1}.txt".format(args.outfile, endosperm)
		write_dict_to_file (all_sites, outfilename)
		
		# annotate
		endo_sites_annotated = get_annotation (set(all_sites.keys()), parent_sites_annotated)
		
		# sum up the counts per gene and write to file
		pergene_counts = count_per_gene (all_sites, endo_sites_annotated)
		outfilename = "{0}.endosperm_{1}_counts_per_gene.txt".format(args.outfile, endosperm)
		write_genecounts_to_file (pergene_counts, outfilename)
		
#	write_output (both_hom, rc_header, "{0}.both_plants_hom_for_diff_alleles.txt".format(args.outfile))
#	write_output (plant_1_heto, rc_header, "{0}.plant_1_heto.txt".format(args.outfile))
#	write_output (plant_2_heto, rc_header, "{0}.plant_2_heto.txt".format(args.outfile))
	
	print "Done!"
	
	
if __name__ == "__main__": 
	main(sys.argv[1:])
