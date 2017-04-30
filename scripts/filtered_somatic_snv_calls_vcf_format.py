import sys
import os
import re

# Input file
infile = sys.argv[1]

# Output file
outfile = open(sys.argv[2], 'w')


# Reference files
birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
header_file = "./data/somaticseq_vcf/738-1_S1_somaticseq_snv_vep.vcf"

# Create an array of vcf files
vcf_file = {}
# Iterate over samples labels (birds) and create array
for bird in open(birds_file):
	bird = bird.rstrip()
	#vcf_file.append(bird)
	vcf_file[bird] = "./data/somaticseq_vcf/" + bird + "_somaticseq_snv_vep.vcf"

# Create a header for the vep vcf file
# Note: I checked all headers from tumor vcf files, they are all the same, can choose one arbitrarily
for header_line in open(header_file):
	if header_line[0] == '#':
		outfile.write(header_line)

# Start at each line in the final filtered variants file and pull appropriate calls from annotated vcf file
for in_line in open(infile):
	if in_line[0] != '#':
		in_line = in_line.rstrip()
		in_cols = in_line.split('\t')
		in_chr = in_cols[0]
		in_pos = in_cols[1]
		in_ref = in_cols[2]
		in_alt = in_cols[3]
		in_snv = in_chr + in_pos + in_ref + in_alt
		in_sample = in_cols[9]
		in_sample_num = in_sample.count('|') + 1
		#print('\n' + in_line)
		#print(in_sample)
		#print(in_sample_num)
		# Note: Range treats numbers as zero-based
		sample = []
		for n in range(in_sample_num):
			print(n)
			print(in_sample)
			print(in_sample.split('|')[n])
			sample.append(n)
			sample[n] = in_sample.split('|')[n]
			print('sample[' + str(n) + ']: ' + sample[n])
			vcf_file2read = vcf_file[sample[n]]
			#print(vcf_file2read)
			for vcf_line in open(vcf_file2read):
				#vcf_line = vcf_line.rstrip()
				if vcf_line[0] != '#': 
					vcf_cols = vcf_line.split('\t')
					vcf_chr = vcf_cols[0]
					vcf_pos = vcf_cols[1]
					vcf_ref = vcf_cols[3]
					vcf_alt = vcf_cols[4]
					vcf_snv = vcf_chr + vcf_pos + vcf_ref + vcf_alt
					if in_snv == vcf_snv:
						outfile.write(vcf_line)
