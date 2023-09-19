import argparse

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--my_samples')
	parser.add_argument('--all_samples')
	parser.add_argument('--out_samples')

	return parser.parse_args()

def main():
	
	## remove bad ids from sample file before generating vcf
	## bad ids = duplicate ids, and ids not in chi file
	
	args = get_args()
	samples_f = args.my_samples
	avail_f = args.all_samples
	out_f = args.out_samples

	my_samples = make_my_samples(samples_f)

	avail_samples = check_available_samples(my_samples, available_f)

	write_avail_samples(avail_samples, out_f)
	
def make_my_samples(samples_file):
	my_samples = []

	s_open = open(samples_file, 'r')
	for line in s_open:	
		L = line.strip()
		if L not in my_samples:
			my_samples.append(L)
	s_open.close()
	reutrn my_samples

def check_available_samples(my_samples, available_file):
	avail_samples = []
	
	a_open = open(available_file, 'r')
	for line in a_open:
		L = line.strip()
		if L in my_samples:
			if L not in avail_samples:
				avail_samples.append(L)
	a_open.close()
	return avail_samples

def write_avail_samples(avail_samples, out_file)
	
	o_open = open(out_file, 'w')
	for avs in avail_samples:
		o_open.write(avs)
		o_open.write('\n')
	o_open.close()


if __name__ == '__main__':
	main()
