import os
import sys
import argparse
import numpy as np
import faiss

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--emb_dir', type=str, default='/Users/kristen/PycharmProjects/faiss/embeddings/')
	parser.add_argument('--idx_dir', type=str, default='/Users/kristen/PycharmProjects/faiss/indexes/')
	parser.add_argument('--k', type=int, default=20)
	parser.add_argument('--out_dir', type=str, default='/Users/kristen/PycharmProjects/faiss/results/')
	parser.add_argument('--query_samples', type=str, default='/Users/kristen/PycharmProjects/faiss/query_samples.txt')
	return parser.parse_args()

def main():
	args = parse_args()
	emb_dir = args.emb_dir
	idx_dir = args.idx_dir
	out_dir = args.out_dir
	k = args.k
	query_samples = args.query_samples
	
	# read query samples
	query_samples_list = read_query_samples(query_samples)
	
	# search all segments for each query sample
	for seg_idx in os.listdir(idx_dir):
		segment = seg_idx.split('.')[1]
		# get embeddings for segment
		embedding_file = emb_dir + 'segment.' + segment + '.txt'
		segment_embeddings_dict = read_embeddings(embedding_file)
		# open results file and clear contents
		results_file = out_dir + 'segment.' + segment + '.results.txt'
		open(results_file, 'w').close()
		# search faiss index for each query sample
		for query_sample in query_samples_list:
			# get query sample embedding
			query_sample_embedding = segment_embeddings_dict[query_sample]
			# search faiss index
			D, I = search_faiss_index(idx_dir + seg_idx, query_sample_embedding, k)
			# write results to file
			write_results(results_file, query_sample, segment_embeddings_dict, D, I)

def search_faiss_index(index_file, query_sample_embedding, k):
	# read faiss index
	index = faiss.read_index(index_file)
	# search faiss index
	D, I = index.search(np.array([query_sample_embedding]), k)
	return D, I

def write_results(results_file, query_sample, segment_embeddings_dict, D, I):
	# write faiss search results to file
	o = open(results_file, 'a')
	o.write('Query: ' + query_sample + '\n')
	# write results
	# sample_hap, distance, index
	for i in range(len(I[0])):
		sample_hap = list(segment_embeddings_dict.keys())[I[0][i]]
		distance = D[0][i]
		index = I[0][i]
		o.write(sample_hap + '\t' + str(distance) + '\n')
	o.write('\n')
	o.close()

def read_embeddings(emb_file):
	# return array of embeddings
	embeddings = {}
	
	# open embedding file
	for line in open(emb_file, 'r'):
	    # split line
	    line = line.split(' ')
	    # get segment name
	    sample_hap = line[0]
	    # get segment embedding
	    segment_embedding = np.array([float(i) for i in line[1:]])
	    # append to embeddings
	    embeddings[sample_hap] = segment_embedding
	return embeddings

def read_query_samples(query_samples):
	# return list of query samples
	query_samples_list = []
	for line in open(query_samples, 'r'):
		query_samples_list.append(line.strip())
	return query_samples_list

if __name__ == '__main__':
	main()
