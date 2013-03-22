#!/usr/bin/env python

import argparse
import os
import sys
import random
import copy
from operator import mul


def genModelInteraction(values, indices):
	return reduce(mul, (values[i] for i in indices))
	
def modelInteraction(indices):
	return lambda x: genModelInteraction(x, indices)

def getModelResult(line, model):
	"""
	Returns the model result for a line
	"""
	return model([float(l) for l in line.split()])

def printOneFile(f_in, model, name_idx, header=False):
	"""
	Prints the output in one file mode
	"""
	fn_in = os.path.splitext(f_in.name)
	f_out = file(fn_in[0] + fn_in[1] + ".model", 'w')
	
	if header:
		expr_vars = [(v,k) for (k,v) in name_idx.iteritems() if v >=0]
		expr_vars.sort()
		expr_str = ' '.join([t[1] for t in expr_vars])
		snp_vars = [(v,k) for (k,v) in name_idx.iteritems() if v < 0]
		snp_vars.sort()
		snp_str = ' '.join([t[1] for t in snp_vars])
		
		l = "#"
		while l.startswith("#"):
			l = f_in.next().strip()

		mid_names = []				
		if len(l.strip().split()) > len(expr_vars) + len(snp_vars):
			mid_names = [str(i) for i in xrange(len(expr_vars)+1, len(l.strip()) - len(snp_vars) + 1)]
					
		print >> f_out, "#result",
		print >> f_out, expr_str,
		if mid_names:
			print >> f_out, ' '.join(mid_names),
		print >> f_out, snp_str
	
		print >> f_out, getModelResult(l, model),
		print >> f_out, l
	
	for l in f_in:
		l = l.strip()
		print >> f_out, getModelResult(l, model),
		print >> f_out, l
		
	f_out.close()
	
def printTwoFile(f_in, model, name_idx, header=False):
	"""
	Prints the output in two file mode
	"""
	fn_in = os.path.splitext(f_in.name)
	f_out_model = file(fn_in[0] + fn_in[1] + ".model", 'w')
	f_out_expr = file(fn_in[0] + fn_in[1] + ".expr", 'w')
	
	expr_vars = [(v,k) for (k,v) in name_idx.iteritems() if v >=0]
	expr_vars.sort()
	expr_str = ' '.join([t[1] for t in expr_vars])
	snp_vars = [(v,k) for (k,v) in name_idx.iteritems() if v < 0]
	snp_vars.sort()
	snp_str = ' '.join([t[1] for t in snp_vars])
	
	if header:	
		l = "#"
		while l.startswith("#"):
			l = f_in.next().strip()
		
		mid_names = []
		if len(l.split()) > len(expr_vars) + len(snp_vars):
			mid_names = [str(i) for i in xrange(len(expr_vars)+1, len(l.split()) - len(snp_vars) + 1)]
		
		print >> f_out_model, "#result",
		if mid_names:
			print >> f_out_expr, "#" + expr_str,
			print >> f_out_expr, ' '.join(mid_names)
		else:
			print >> f_out_expr, "#" + expr_str
		
		print >> f_out_model, snp_str
	
		print >> f_out_model, getModelResult(l, model),
		print >> f_out_model, l.split()[len(snp_vars):]
		print >> f_out_expr, l.split()[:len(snp_vars)]
	
	for l in f_in:
		l = l.strip()
		print >> f_out_model, getModelResult(l, model),
		print >> f_out_model, ' '.join(l.split()[len(snp_vars):])
		print >> f_out_expr, ' '.join(l.split()[:len(snp_vars)])
		
	f_out_model.close()
	f_out_expr.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Calculate a model output given a model definition")
	
	parser.add_argument("--result", "-r", nargs="*", type=argparse.FileType('r'), help="Result File(s)")
	parser.add_argument("--map", "-m", help="Filename of names of SNP variables", type=argparse.FileType('r'))
	parser.add_argument("--expression", "-x", help="Filename of names of expression variables", type=argparse.FileType('r'))
	parser.add_argument("--model", "-d", help="Filename of the model specification", required=True, type=argparse.FileType('r'))
	parser.add_argument("--format", "-f", help="Format of the output", choices=["one_file", "two_file"], default="one_file")
	parser.add_argument("--header", "-a", help="Include header in output", action="store_true")
	parser.add_argument("--stddev", "-s", type=float, help="Standard deviation", default=1)
	parser.add_argument("--seed", "-e", type=int, help="Random number seed", default=None)
	
	options = parser.parse_args()
	
	# Set the random number seed
	if options.seed is not None:
		random.seed(options.seed)
	
	name_idx = {}
	# Read the map and expression files and convert them to indices
	# Note: map count from the end, expression count from the beginning
	snp_list = []
	if options.map:
		for l in options.map:
			l = l.strip()
			if not l.startswith("#"):
				snp_list.append(l.split()[1])
		
		name_idx = dict(zip(snp_list,(i-len(snp_list) for i in xrange(len(snp_list)))))
		
	i = 0
	if options.expression:
		for l in options.expression:
			l = l.strip()
			if not l.startswith("#"):
				name_idx[l] = i
				i += 1
	
	model_fns = []
	# Read the model definition file
	# Note: if not found in the index above, try converting to an integral index
	for l in options.model:
		l = l.strip()
		if not l.startswith("#"):
			(coeff, models) = l.split('\t')
			models = models.split(',')
			coeff = float(coeff)
			idx_list = []
			for m in models:
				m = m.strip()
				if m in name_idx:
					idx_list.append(name_idx[m])
				else:
					idx_list.append(int(m)-1*(int(m)>0))
					
			model_fns.append(modelInteraction(idx_list))
			
	final_model = lambda x: random.gauss(sum((f(x) for f in model_fns)), options.stddev)
	
	for r in options.result:
		if options.format=="one_file":
			printOneFile(r, final_model, name_idx, options.header)
		elif options.format == "two_file":
			printTwoFile(r, final_model, name_idx, options.header)
