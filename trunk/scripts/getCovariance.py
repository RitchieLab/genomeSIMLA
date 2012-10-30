#!/usr/bin/env python

import os
import sys
import argparse
import random
from numpy import cov

def readGEOFile(f, N, inc_set, sel_method):
	"""
	This method reads a gene omnibus file and returns the data
	
	f: the GEO file
	N: # of lines to select
	inc_set: set of ids to include by force
	sel_method: function that dictates the selection method
	"""
	
	data_keys = set([])
	for l in f:
		l = l.strip()
		if l.startswith("#"):
			l = l.strip("#")
			if l.split("=")[1].strip().startswith("Value"):
				data_keys.add(l.split("=")[0].strip())
		
		if l == "!dataset_table_begin":
			break
	
	l = f.next()
	curr_idx = 0
	id_idx = -1;
	data_idx = []
	for d in l.split():
		if d in data_keys:
			data_idx.append(curr_idx)
		elif d == "ID_REF":
			id_idx = curr_idx
		
		curr_idx += 1 
	
	line_list = sel_method(f, N, "\t", inc_set, id_idx, "!")
	data_list = [[]]*len(line_list)
	
	if id_idx == -1:
		raise Exception("Could not find identifier column")
	
	if len(data_idx) == 0:
		raise Exception("Could not find data columns")
	
	curr_idx = 0
	id_list = []
	for l in line_list:
		if not l:
			raise Exception("Required Line not found or not enough lines")
		
		data_str = l.split("\t")
		
		id_list.append(data_str[id_idx])
		missing_idx = []
		n_samples = 0
		total = 0
		data_list[curr_idx] = [0]*len(data_idx)
		curr_data_idx = 0
		for d_idx in data_idx:
			try:
				data_list[curr_idx][curr_data_idx] = float(data_str[d_idx])
				n_samples += 1
				total += data_list[curr_idx][curr_data_idx]
			except ValueError:
				missing_idx.append(curr_data_idx)
			
			curr_data_idx += 1
		
		if n_samples == 0:
			raise Exception("Sample found with no data")
		
		avg = total / n_samples
		
		for m_idx in missing_idx:
			data_list[curr_idx][m_idx] = avg
			
		curr_idx += 1
		
	return (id_list, data_list)
	
def selectRandom(f, N, sep, inc_set, id_idx, comment="#"):
	"""
	This method selects N random lines from the file
	
	f is the file object
	N is the total number of lines to select
	sep is the field separator of the line
	inc_set is the set of ids to include
	id_idx is the index of the ids to test for inclusion
	"""
	r = random.Random()
	ret_list = [''] * N
	
	n_lines = 0
	n_choice = N - len(inc_set)
	choice_idx = n_choice
	for l in f:
		l = l.strip()
		if not l.startswith(comment):
			id = l.split(sep)[id_idx]
			# If this is one of the mandatory selections, put it toward the end
			if id in inc_set:
				ret_list[choice_idx] = l
				choice_idx += 1
			# Otherwise, select randomly (also prevents duplicates)
			else:
				idx = r.randint(0,n_lines)
				if idx < n_choice:
					# Fill the reservoir first if not full
					if n_lines < n_choice:
						ret_list[n_lines] = l
					# Then, we have to replenish individual bins
					else:
						ret_list[idx] = l
				n_lines += 1
			
	return ret_list
	
def selectFirst(f, N, sep, inc_set, id_idx, comment="#"):
	"""
	This method selects the first N lines from the file
	"""
	ret_list = ['']*N
	idx = 0
	n_choice = N - len(inc_set)
	choice_idx = n_choice
	for l in f:
		l = l.strip()
		if not l.startswith(comment):
			id = l.split(sep)[id_idx]
			# If this is one of the mandatory selections, put it toward the end
			if id in inc_set:
				ret_list[choice_idx] = l
				choice_idx += 1
			# Otherwise, select randomly (also prevents duplicates)
			elif idx < n_choice:
				ret_list[idx] = l
				idx += 1
		
	return ret_list
	
def selectLast(f, N, sep, inc_set, id_idx, comment="#"):
	"""
	This method selects the last N lines from the file
	"""
	ret_list = ['']*N
	idx = 0
	n_choice = N - len(inc_set)
	choice_idx = n_choice
	for l in f:
		l = l.strip()
		if not l.startswith(comment):
			id = l.split(sep)[id_idx]
			# If this is one of the mandatory selections, put it toward the end
			if id in inc_set:
				ret_list[choice_idx] = l
				choice_idx += 1
			# Otherwise, select into a ring buffer
			else:
				ret_list[idx] = l
				idx = (idx + 1) % n_choice

	return ret_list

def calcCovariance(data_in):
	"""
	Get the covariance from a list of lists of data
	"""
	return cov(data_in)

if __name__ == "__main__":
	
	reader_dict = {"SOFT": readGEOFile}
	select_dict = {"random": selectRandom,
	               "first" : selectFirst,
	               "last" : selectLast}
	
	# define the options
	parser = argparse.ArgumentParser(description="Calculate a covariance matrix from a given gene expression file.\n\n" + 
	"The covariance matrix is printed on stdout, and the IDs of the data used to generate the covariance matrix are printed on stderr.")
	
	# name of the data file
	parser.add_argument("--file","-f", nargs="?", 
		type=argparse.FileType('r'), default=sys.stdin,
		help="Filename to read gene expression data from (if not given, read from standard input)")
	
	# format of the data file
	parser.add_argument("--type","-t", choices=reader_dict.keys(), default="SOFT", help="Format of gene expression data")
	
	# selection strategy
	parser.add_argument("--selection","-s", choices=select_dict.keys(), default="random", help="Method of selection")
	
	# number to be selected
	parser.add_argument("--number","-n", type=int, default=50, help="Number of gene expression variables to calculate covariance for")
	
	# Any that are required to be selected	
	parser.add_argument("--required","-r", nargs="*", help="List of gene expression variables to include in the analysis", default=[])	
	
	options = parser.parse_args()
	# get the type of the file
	read_fn = reader_dict[options.type]
		
	# get the name of the file
	
	# get the style of selection
	sel_fn = select_dict[options.selection]
	
	# use the correct reading function to get the data
	(cov_names, cov_data) = read_fn(options.file, options.number, set(options.required), sel_fn)
	
	# Print the covariance names:
	print >> sys.stderr, "# IDs included in covariance analysis:"
	for n in cov_names:
		print >> sys.stderr, n
	
	# Now, calculate and print the covariance
	c = calcCovariance(cov_data)
	for l in c:
		print '\t'.join(("%f" % i for i in l))
		
	 
	