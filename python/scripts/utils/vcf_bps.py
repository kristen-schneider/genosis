import sys, os, gzip

vcf_file = gzip.open(sys.argv[1]) #.bed
outfile = open(sys.argv[2], "w") #output style: [rs] [pos] [genetic pos]

vcf_pos = []
vcf_chr = []
vcf_info = []

# make list of all bp in vcf
print('reading vcf...')
for line in vcf_file:
    #ignore header
    if line.startswith(b'#'):
        continue
    else:
        line = line.strip().split()
        chrm = int(line[0].replace(b'chr',b''))
        pos = int(line[1])
        info = line[2].decode()
        
        line = str(chrm) + '\t' + str(pos) + '\t' + str(info) + '\n'
        outfile.write(line)
        #rs = line[3]
        #vcf_chr.append(chrm)
        #vcf_pos.append(pos)
        #vcf_info.append(info)

#print('writing out...')
#for vcf_idx in range(len(vcf_pos)):
#    line = str(vcf_chr[vcf_idx] + '\t' + vcf_info[vcf_idx] + '\t' + vcf_pos[vcf_idx] + '\n')
#    outfile_write(line)

