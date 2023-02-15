#:https://github.com/joepickrell/1000-genomes-genetic-maps/blob/master/scripts/interpolate_maps.py
#!/usr/bin/python


import sys, os, gzip

vcf_bps = open(sys.argv[1]) #.bed
map_file = open(sys.argv[2]) #input map file, either the HapMap map or the 1000 Genomes OMNI map
outfile = open(sys.argv[3], "w") #output style: [rs] [pos] [genetic pos]

vcf_chr = []
vcf_pos = []
vcf_info = []
map_bp = []
map_cm = []

# make list of all bp in vcf
print('reading vcf bps...')
for line in vcf_bps:
    line = line.strip().split()
    chrm = int(line[0])
    pos = int(line[1])
    info = line[2]
    vcf_chr.append(chrm)
    vcf_pos.append(pos)
    vcf_info.append(info)

print('reading map file...')
# make list of all bp and cm in map file
for line in map_file:
	line = line.strip().split()
	bp_pos = int(line[3])
	cm_pos = float(line[2])
    #pos = int(line[1]) #uncomment for hapmap input
	#gpos = float(line[2])
	#gpos = float(line[3]) #uncomment for hapmap  input
	map_bp.append(bp_pos)
	map_cm.append(cm_pos)

vcf_idx = 0
map_idx = 0
# iterator over all bp pos in vcf
while vcf_idx < len(vcf_pos):
    
    pos = vcf_pos[vcf_idx]
    #rs = rsin[vcf_idx]
    
    if (pos == map_bp[map_idx]):
		#the 1000 Genomes site was genotyped as part of the map
        out_line = str(vcf_chr[vcf_idx]) + '\t' + str(vcf_info[vcf_idx]) + '\t' + str(map_cm[map_idx]) + '\t' + str(vcf_pos[vcf_idx]) + '\n'
        outfile.write(out_line)
        #print(out_line)
        vcf_idx += 1
    elif(pos < map_bp[map_idx]):
        #current position in interpolation before marker
        if (map_idx == 0):
            #before the first site in the map (genetic position = 0)
            out_line = str(vcf_chr[vcf_idx]) + '\t' + str(vcf_info[vcf_idx]) + '\t' + str(map_cm[map_idx]) + '\t' + str(vcf_pos[vcf_idx]) + '\n'
            outfile.write(out_line)
            #print(out_line)
            #outfile.write(rs + "\t" + pos + "\t" + map_cm[map_idx])
            vcf_idx += 1
        else:
            #interpolate
            prev_map_cm = map_cm[map_idx-1]
            prev_map_pos = map_bp[map_idx]
            
            frac = (float(pos)-float(map_bp[map_idx-1]))/ (float(map_bp[map_idx]) - float(map_bp[map_idx-1]))
            temp_cm = prev_map_cm + frac * (map_cm[map_idx]-prev_map_cm)
            out_line = str(vcf_chr[vcf_idx]) + '\t' + str(vcf_info[vcf_idx]) + '\t' + str(temp_cm) + '\t' + str(vcf_pos[vcf_idx]) + '\n'
            outfile.write(out_line)
            #print(out_line)
            vcf_idx += 1
    
    elif (pos > map_bp[map_idx]):
        #current position in interpolation after marker
        if (map_idx == len(map_bp)-1):
            #after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
            out_line = str(vcf_chr[vcf_idx]) + '\t' + str(vcf_info[vcf_idx]) + '\t' + str(map_cm[map_idx]) + '\t' + str(vcf_pos[vcf_idx]) + '\n'
            #out_line = vcf_chr[vcf_idx] + '\t' + vcf_info[vcf_idx] + '\t' + map_cm[map_idx] + '\t' + vcf_pos[vcf_idx]
            outfile.write(out_line)
            #print(out_line)
            vcf_idx += 1
        else:
           #increment the marker
           map_idx += 1
