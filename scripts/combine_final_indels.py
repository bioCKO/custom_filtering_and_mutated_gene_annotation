import sys
import os
import re

# Output file
outfile = open("./results/somatics_indels_final.txt", 'w')

# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"
final_somatic_indel_files = '/home/proj/MDW_genomics/steepale/pathway_analysis/data/final_somatic_indel_files.txt'

# Write header to outfile
outfile.write('#CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'MUT' + '\t' + 'IMPACT' + '\t' + 'SYMBOL' + '\t' + 'GENE_ID' + '\t' + 'TSN' + '\t' + 'SAMPLE' + '\t' + 'VAC' + '\t' + 'VAF' + '\n')

# Iterate over samples labels (birds) and create dictionary
final_indel_file = {}
for bird in open(tumor_birds_file):
	bird = bird.rstrip()
	#vcf_file.append(bird)
	final_indel_file[bird] = "./data/all_nonsyn_indels_" + bird + "_final.txt"

# Define empty variables
#final_output = None

other_file_list = ['./data/all_nonsyn_indels_738-1_S1_final.txt',
'./data/all_nonsyn_indels_741-1_S2_final.txt',
'./data/all_nonsyn_indels_756-3_S3_final.txt',
'./data/all_nonsyn_indels_766-1_S4_final.txt',
'./data/all_nonsyn_indels_777-3_S14_final.txt',
'./data/all_nonsyn_indels_787-2_S15_final.txt',
'./data/all_nonsyn_indels_788-1_S16_final.txt',
'./data/all_nonsyn_indels_794-1_S17_final.txt',
'./data/all_nonsyn_indels_798-1_S5_final.txt',
'./data/all_nonsyn_indels_833-1_S6_final.txt',
'./data/all_nonsyn_indels_834-2_2_S12_final.txt',
'./data/all_nonsyn_indels_834-2_S7_final.txt',
'./data/all_nonsyn_indels_835-1_S18_final.txt',
'./data/all_nonsyn_indels_841-3_S19_final.txt',
'./data/all_nonsyn_indels_842-2_2_S25_final.txt',
'./data/all_nonsyn_indels_842-2_S20_final.txt',
'./data/all_nonsyn_indels_855-1_S8_final.txt',
'./data/all_nonsyn_indels_863-1_S9_final.txt',
'./data/all_nonsyn_indels_884-2_S21_final.txt',
'./data/all_nonsyn_indels_901-2_2_S26_final.txt',
'./data/all_nonsyn_indels_901-2_S22_final.txt',
'./data/all_nonsyn_indels_906-1_S23_final.txt',
'./data/all_nonsyn_indels_911-1_2_S13_final.txt',
'./data/all_nonsyn_indels_911-1_S24_final.txt',
'./data/all_nonsyn_indels_918-3_S10_final.txt',
'./data/all_nonsyn_indels_927-2_S11_final.txt']

# Grab unique values and combine variable
# Start going through all the samples one by one
for bird in open(tumor_birds_file):
	bird = bird.rstrip()
	# For each line in final SNVs file
	for indel_line in open(final_indel_file[bird]):
		indel_line = indel_line.rstrip()
		if indel_line[0] != '#':
			indel_col = indel_line.split('\t')
			indel_chr = indel_col[0]
			final_chr = indel_chr
			indel_pos = indel_col[1]
			final_pos = indel_pos
			indel_ref = indel_col[2]
			final_ref = indel_ref
			indel_alt = indel_col[3]
			final_alt = indel_alt
			indel_mut = indel_col[4]
			final_mut = indel_mut
			indel_impact = indel_col[5]
			final_impact = indel_impact
			indel_symbol = indel_col[6]
			final_symbol = indel_symbol
			indel_geneid = indel_col[7]
			final_geneid = indel_geneid
			indel_tsn = int(indel_col[8])
			final_tsn = indel_tsn
			indel_sample = indel_col[9]
			final_sample = indel_sample
			indel_vac = indel_col[10]
			final_vac = indel_vac
			indel_vaf = indel_col[11]
			final_vaf = indel_vaf
			#final_output = None
			#print(final_tsn)
			for n in range(26):
				for other_line in open(other_file_list[n]):
					other_line = other_line.rstrip()
					#print('other_line: ' + other_line)
					if other_line[0] != '#':
						other_col = other_line.split('\t')
						other_chr = other_col[0]
						other_pos = other_col[1]
						other_ref = other_col[2]
						other_alt = other_col[3]
						other_mut = other_col[4]
						other_impact = other_col[5]
						other_symbol = other_col[6]
						other_geneid = other_col[7]
						other_tsn = int(other_col[8])
						other_sample = other_col[9]
						other_vac = other_col[10]
						other_vaf = other_col[11]
						if indel_chr == other_chr and indel_pos == other_pos and indel_alt == other_alt and indel_sample != other_sample:
							final_output = 'Won'
							#final_chr = indel_chr
							#final_pos = indel_pos
							#final_ref = indel_ref
							#final_alt = indel_alt
							#final_mut = indel_mut
							#final_impact = indel_impact
							#final_symbol = indel_symbol
							#final_geneid = indel_geneid
							final_tsn = final_tsn + other_tsn
							final_sample = final_sample + '|' + other_sample 
							final_vac = final_vac + '|' + other_vac
							final_vaf = final_vaf + '|' + other_vaf
			#if final_output != None:
			outfile.write(final_chr + '\t' + final_pos + '\t' + final_ref + '\t' + final_alt + '\t' + final_mut + '\t' + final_impact + '\t' + final_symbol + '\t' + final_geneid + '\t' + str(final_tsn) + '\t' + final_sample + '\t' + final_vac + '\t' + final_vaf + '\n')

