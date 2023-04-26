import os
import argparse
import numpy as np
import faiss

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--emb_dir', type=str)
	parser.add_argument('--idx_dir', type=str)
	return parser.parse_args()

def main():
	args = parse_args()
	emb_dir = args.emb_dir
	idx_dir = args.idx_dir
	
	# for a directory of gt indexes
	for gt_encoding in os.listdir(emb_dir):
		if gt_encoding.endswith('.txt'):
			print('Encoding file: {}'.format(gt_encoding))
			
			# read data from file
			embeddings_numpy = read_embeddings(emb_dir + gt_encoding)
			
			# build faiss index for l2 distance and write to file
			l2_index = build_l2_index(embeddings_numpy)
			base_name = '.'.join(gt_encoding.split('.')[0:2])
			l2_index_file = idx_dir + base_name + '.index.l2'
			faiss.write_index(l2_index, l2_index_file)
		break

def read_embeddings(emb_file):
	# return array of embeddings
	embeddings = []
	
	# open embedding file
	for line in open(emb_file, 'r'):
		# split line
		line = line.split(' ')
		# get segment name
		sample_hap = line[0]
		# get segment embedding
		segment_embedding = np.array([float(i) for i in line[1:]])
		# append to embeddings
		embeddings.append((sample_hap, segment_embedding))

	# convert data to numpy array and reshape
	data_array = np.array([i[1] for i in embeddings])
	data_array = data_array.reshape(data_array.shape[0], data_array.shape[1])
	return data_array

def build_l2_index(data_array):
	# build faiss index for l2 distance
	index = faiss.IndexFlatL2(data_array.shape[1])
	index.add(data_array)
	return index

if __name__ == '__main__':
	main()
