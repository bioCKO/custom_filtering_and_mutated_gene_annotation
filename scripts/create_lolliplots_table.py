import sys
import os
import re

# General Strategy of script:
# Take the final somatic snvs and indels in vep annotation format and create a genvisr table for lolliplots

# Input somatic SNV file
#snv_infile = open(sys.argv[1], 'r')
# Input somatic indel file 
#indel_infile = open(sys.argv[2], 'r')

# Reference files
# final somatic SNVs and Indels
snvs_indels_sample_file = "./results/somatic_snvs_and_indels_final.txt"

# read the combined somatic snv and indel file
muts_file = "./results/somatic_snvs_indels_combined_vep_no_header.vcf"
# Output file
outfile = open("./results/somatic_snvs_indels_genvisr.int", 'w')

# Write header to outfile
# Note: I have to have certain headers which is really annoying
outfile.write('gene' + '\t' + 'amino_acid_change' + '\t' + 'transcript_name' + '\t' + 'protein_id' + '\t' + 'Impact' + '\t' + 'Sample' + '\n')

# Loop through the snv and indel final muts file and grab samples
for snv_indel_line in open(snvs_indels_sample_file):
	if snv_indel_line[0] != '#':
		snv_indel_line = snv_indel_line.rstrip()
		snv_indel_col = snv_indel_line.split('\t')
		snv_indel_chr = snv_indel_col[0]
		snv_indel_pos = snv_indel_col[1]
		snv_indel_ref = snv_indel_col[2]
		snv_indel_alt = snv_indel_col[3]
		snv_indel_symbol = snv_indel_col[6]
		snv_indel_gene_id = snv_indel_col[7]
		#print('file 1 line')
		snv_indel_mutation = snv_indel_chr + '\t' + snv_indel_pos + '\t' + snv_indel_ref + '\t' + snv_indel_alt
		snv_indel_samples = snv_indel_col[9]
		for muts_line in open(muts_file):
			if muts_line[0] != '#':
				#muts_line = muts_line.rstrip()
				#print(muts_line)
				muts_cols = muts_line.split('\t')
				muts_chr = muts_cols[0]
				muts_pos = muts_cols[1]
				muts_ref = muts_cols[3]
				muts_alt = muts_cols[4]
				#print('file 2 line')
				muts_mutation = muts_chr + '\t' + muts_pos + '\t' + muts_ref + '\t' + muts_alt
				#if snv_indel_chr == muts_chr:
					#print('connection')
				muts_info = muts_cols[7]
				muts_format = muts_cols[8]
				muts_normal = muts_cols[9]
				muts_tumor = muts_cols[10]
				#print(vep_snv)
				#print(vep_sample)
				if re.search('SOMATIC', muts_info.split(';')[0]):
					muts_somat = muts_info.split(';')[0]
				else:
					muts_somat = 'NA'
				if re.search('MVJSDU', muts_info.split(';')[0]):
					muts_tools = muts_info.split(';')[0]
				elif re.search('MVJSDU', muts_info.split(';')[1]):
					muts_tools = muts_info.split(';')[1]
				if re.search('NUM_TOOLS', muts_info.split(';')[1]):
					muts_tool_num = muts_info.split(';')[1]
				elif re.search('NUM_TOOLS', muts_info.split(';')[2]):
					muts_tool_num = muts_info.split(';')[2]
				if re.search('CSQ=', muts_info.split(';')[2]):
					muts_info_ann = muts_info.split(';')[2].split('SQ=')[1]
					muts_info_ann_num = muts_info_ann.count(',') + 1
					info_ann = []
				elif re.search('CSQ=', muts_info.split(';')[3]):
					muts_info_ann = muts_info.split(';')[3].split('SQ=')[1]
					muts_info_ann_num = muts_info_ann.count(',') + 1
					info_ann = []
				for n in range(muts_info_ann_num):
					info_ann.append(n)
					info_ann[n] = muts_info_ann.split(',')[n]
					info_ann2read = info_ann[n]
					info_cols = info_ann2read.split('|')
					muts_allele = info_cols[0]
					muts_cons = info_cols[1]
					muts_impact = info_cols[2]
					muts_symbol = info_cols[3]
					muts_geneid = info_cols[4]
					muts_feat_type = info_cols[5]
					muts_feature = info_cols[6]
					muts_biotype = info_cols[7]
					muts_exon = info_cols[8]
					muts_intron = info_cols[9]
					muts_HGVSc = info_cols[10]
					muts_HGVSp = info_cols[11]
					muts_cDNA_pos = info_cols[12]
					muts_CDS_pos = info_cols[13]
					muts_protein_pos = info_cols[14]
					muts_aminos = info_cols[15]
					muts_codons = info_cols[16]
					muts_existing_var = info_cols[17]
					muts_distance = info_cols[18]
					muts_strand = info_cols[19]
					muts_flags = info_cols[20]
					muts_symbol_source = info_cols[21]
					muts_HGNC_ID = info_cols[22]
					muts_tsl = info_cols[23]
					muts_appris = info_cols[24]
					muts_ccds = info_cols[25]
					muts_ensp = info_cols[26]
					muts_swissprot = info_cols[27]
					muts_trembl = info_cols[28]
					muts_uniparc = info_cols[29]
					muts_sift = info_cols[30]
					muts_domains = info_cols[31]
					muts_hgvs_offset = info_cols[32]
					if muts_symbol == '':
						muts_symbol = 'NA'
					if snv_indel_symbol == 'NA':
						output_symbol = snv_indel_gene_id
					else:
						output_symbol = snv_indel_symbol
					#print(muts_HGVSp)
					if muts_HGVSp == '':
						output_HGVSp = 'remove'
					else:
						half2_HGVSp = muts_HGVSp.split(':')[1]
						output_HGVSp = half2_HGVSp
						#p_HGVSp = half2_HGVSp[0:2]
						#AA1_HGVSp = half2_HGVSp.split('p.')[1][0:3]
						#if re.findall(r'3D', half2_HGVSp):
						#	last_part_HGVSp = half2_HGVSp[-3:]
						#	pos_HGVSp = re.findall(r'\d+', half2_HGVSp[:-3])
						#elif:
						#else:
						#	pos_HGVSp = re.findall(r'\d+', half2_HGVSp)
						#output_HGVSp = muts_HGVSp
					if muts_feature == '':
						output_feature = 'NA'
					else:
						output_feature = muts_feature
					if muts_ensp == '':
						output_ensp = 'NA'
					else:
						output_ensp = muts_ensp
					if muts_sift == '':
						output_sift = 'Not_App'
					else:
						output_sift = muts_sift
					output_samples = snv_indel_samples.replace('|', ', ')
					if snv_indel_mutation == muts_mutation:
						outfile.write(output_symbol + '\t' + output_HGVSp + '\t' + output_feature + '\t' + output_ensp + '\t' + output_sift + '\t' + output_samples + '\n')

#Close the output file
outfile.close()
