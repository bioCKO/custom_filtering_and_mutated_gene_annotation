import sys
import os

muts_infile = "./data/all_nonsyn_snvs_741-1_S2_final.txt"

gatk_infile = "./data/GATK_mut_calls/741-1_S2_raw_snps_indels_genotyped.g.vcf"

outfile = open("./data/all_nonsyn_snvs_741-1_S2_final.g.vcf", 'w')

# Write a header from the GATK input file to the output file
for gatk_line in open(gatk_infile):
	if gatk_line[0] == '#':
		outfile.write(gatk_line)

# Search both files for common variants
for muts_line in open(muts_infile):
	muts_line = muts_line.rstrip()
	if muts_line[0] != '#':
		muts_col = muts_line.split('\t')
		muts_chr = muts_col[0]
		muts_pos = muts_col[1]
		muts_ref = muts_col[2]
		muts_alt = muts_col[3]
		for gatk_line in open(gatk_infile):
			gatk_line = gatk_line.rstrip()
			if gatk_line[0] != '#':
				gatk_col = gatk_line.split('\t')
				gatk_chr = gatk_col[0]
				gatk_pos = gatk_col[1]
				gatk_ref = gatk_col[3]
				gatk_alt = gatk_col[4]
				if muts_chr == gatk_chr and muts_pos == gatk_pos and muts_ref == gatk_ref and muts_alt == gatk_alt:
					outfile.write(gatk_line + '\n')
