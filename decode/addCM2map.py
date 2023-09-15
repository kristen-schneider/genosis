import pandas as pd
from scipy.interpolate import interpld
import sys
import os
import argparse

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(--full_map, type=str)
	parser.add_argument(--cmap, type=str)

def main():
	args = parse_args()
	full_map = args.full_map
	cmap = args.cmap
	
	gdf = pd.read_csv(full_map, sep='\t', header=None)
	
	fun - interpld(gdf[2], gdf[6], fill_value='extrapolate')

	mdf = pd.read_csv(cmap, sep='\t')
	mdf['cM'] = fun(mdf.pos)
	
	mdf.to_csv(sys.stdout, index=False, sep='\t')

if __name__ == '__main__':
	main()
