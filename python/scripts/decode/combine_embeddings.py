# sample chr seg embedding
import sys

all_embedding_paths = sys.argv[1]

def main():
	all_embedding_files = read_list_of_files(all_embedding_paths)
	for embedding_file in all_embedding_files:
		write_embeddings(embeddings_all)


def write_embeddings(embeddings_all):
	with open(embeddings_all, 'r') as f:
		for line in f:
			L = line.strip().split()
			sample = L[0]
			chrm = L[1]
			segment = L[2]
			embedding = L[3:]
			# append sample and embedding to segment file
			segment_file = embeddings_dir + 'chrm' + chrm + '.segment' + segment + '.emb'
			with open(segment_file, 'a') as f:
				f.write(sample + ' ' + ' '.join(embedding) + '\n')


def read_list_of_files(all_embedding_paths):
	all_embedding_files = []
	o = open(all_embedding_paths, 'r')
	for line in o:
		L = line.strip().split()
		all_embedding_files.append(L[0])
	o.close()
	return all_embedding_files

main()
