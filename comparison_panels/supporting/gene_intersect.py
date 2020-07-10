#!/usr/bin/env python
import sys

files = sys.argv[1:]
genes = set()

for file in files:
	print file
	with open(file,'r') as f:
		new = set(map(str.strip,f.readlines()))
		if len(genes)==0: 
			genes = new
		else:
			genes = genes.intersection(new)

print genes
print len(genes)
