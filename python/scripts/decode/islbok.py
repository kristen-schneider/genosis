#!/path/to/python2.7

import argparse
import os
import sys
sys.path.insert(1, '/path/to/ibc.py')
import ibc

Q = ibc.IBC(host, port)

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--sample_IDs')
	parser.add_argument('--family_ID')

	return parser.parse_args()

def main():
	args = get_args()
	sample_IDs_pns = args.sample_IDs
	family_ID = args.family_ID
	
	samples = read_samples(sample_IDs_pns)
	for s in samples:
		Q.cmd("father " + s)
		for l in Q.read_lines():
			father = l.strip()
		Q.cmd("mother " + s)
		for l in Q.read_lines():
			mother = l.strip()
		
		print( family_ID + '\t' + s + '\t' + father + '\t' + mother + '\t0\t0' )


def read_samples(sample_IDs_pns):
	samples = []
	
	s_open = open(sample_IDs_pns, 'r')
	for line in s_open:
		L = line.strip().split()
		samples.append(L[0])
	return samples

if __name__ == '__main__':
	main()

