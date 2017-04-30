import sys
import os
import re
import subprocess
from subprocess import check_output
#import logging

# General Strategy of script:
# Script will take a single vep file from one tumor sample and search for each called variant
# across it's own tumor as well as all germline bam files with samtools mpileup. Then filters
# will be applied.

# Filters include:
# Variant must be in tumor bam file from which it came
# Variants must have a variant allele frequency greater than or equal to 0.05
# Tumor file and germline bam pairs much each have a coverage of 4 at variant site
# Variant must not be found in germline files

# Input file
vep_file = open(sys.argv[1], 'r')

# Output file
outfile = open(sys.argv[2], 'w')

# Tumor sample
tbird = sys.argv[3]

# Output to loggin file
#logging.basicConfig(filename='./analysis/vep_filter.log' ,level=logging.DEBUG)

# Reference files
#tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"

# Create a dictionary of germline bam files
germline_bam_file = {}
for gbird in open(germline_birds_file):
	gbird = gbird.rstrip()
	germline_bam_file[gbird] = '/home/proj/MDW_genomics/xu/final_bam/' + gbird + '_Bwa_RG_dedupped_realigned.bam'

# Create a dictionary of tumor bam files
tumor_bam_file = {}
#for tbird in open(tumor_birds_file):
#	tbird = tbird.rstrip()
#tumor_bam_file[tbird] = '/home/proj/MDW_genomics/xu/final_bam/' + tbird + '_Bwa_RG_dedupped_realigned.bam'
tumor_bam_file[tbird] = '/home/proj/MDW_genomics/xu/final_bam/' + tbird + '_Bwa_RG_dedupped_realigned.bam'

# Write header to outfile
outfile.write('#CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'MUT' + '\t' + 'IMPACT' + '\t' + 'SYMBOL' + '\t' + 'GENE_ID' + '\t' + 'TSN' + '\t' + 'SAMPLE' + '\t' + 'VAC' + '\t' + 'VAF' + '\n')


#lin_num = 1
# Create an array of the annotated vcf files
# vcf_file = {}
# Iterate over samples labels (birds) and create array
# for bird in open(tumor_birds_file):
#	bird = bird.rstrip()
	#vcf_file.append(bird)
#	vcf_file[bird] = "./data/somaticseq_vcf/" + bird + "_somaticseq_snv_vep.vcf"
	#Examine each vep file one by one to create a master file
for vep_line in vep_file:
	if vep_line[0] != '#':
		#lin_num = lin_num + 1
		#print(lin_num)
		vep_line = vep_line.rstrip()
		vep_cols = vep_line.split('\t')
		vep_chr = vep_cols[0]
		vep_pos = vep_cols[1]
		vep_ref = vep_cols[3]
		vep_alt = vep_cols[4]
		vep_snv = vep_chr + '\t' + vep_pos + '\t' + vep_ref + '\t' + vep_alt
		vep_info = vep_cols[7]
		vep_format = vep_cols[8]
		vep_normal = vep_cols[9]
		vep_tumor = vep_cols[10]
		vep_sample = tbird.rstrip()
		#print(vep_snv)
		#print(vep_sample)
		if re.search('SOMATIC', vep_info.split(';')[0]):
			vep_somat = vep_info.split(';')[0]
		else:
			vep_somat = 'NA'
		if re.search('MVJSDU', vep_info.split(';')[0]):
			vep_tools = vep_info.split(';')[0]
		elif re.search('MVJSDU', vep_info.split(';')[1]):
			vep_tools = vep_info.split(';')[1]
		if re.search('NUM_TOOLS', vep_info.split(';')[1]):
			vep_tool_num = vep_info.split(';')[1]
		elif re.search('NUM_TOOLS', vep_info.split(';')[2]):
			vep_tool_num = vep_info.split(';')[2]
		if re.search('CSQ=', vep_info.split(';')[2]):
			vep_info_ann = vep_info.split(';')[2].split('SQ=')[1]
			vep_info_ann_num = vep_info_ann.count(',') + 1
			info_ann = []
		elif re.search('CSQ=', vep_info.split(';')[3]):
			vep_info_ann = vep_info.split(';')[3].split('SQ=')[1]
			vep_info_ann_num = vep_info_ann.count(',') + 1
			info_ann = []
		for n in range(vep_info_ann_num):
			info_ann.append(n)
			info_ann[n] = vep_info_ann.split(',')[n]
			info_ann2read = info_ann[n]
			info_cols = info_ann2read.split('|')
			vep_allele = info_cols[0]
			vep_cons = info_cols[1]
			vep_impact = info_cols[2]
			vep_symbol = info_cols[3]
			vep_geneid = info_cols[4]
			vep_feat_type = info_cols[5]
			vep_feature = info_cols[6]
			vep_biotype = info_cols[7]
			vep_exon = info_cols[8]
			vep_intron = info_cols[9]
			vep_HGVSc = info_cols[10]
			vep_HGVSp = info_cols[11]
			vep_cDNA_pos = info_cols[12]
			vep_CDS_pos = info_cols[13]
			vep_protein_pos = info_cols[14]
			vep_aminos = info_cols[15]
			vep_codons = info_cols[16]
			vep_existing_var = info_cols[17]
			vep_distance = info_cols[18]
			vep_strand = info_cols[19]
			vep_flags = info_cols[20]
			vep_symbol_source = info_cols[21]
			vep_HGNC_ID = info_cols[22]
			vep_tsl = info_cols[23]
			vep_appris = info_cols[24]
			vep_ccds = info_cols[25]
			vep_ensp = info_cols[26]
			vep_swissprot = info_cols[27]
			vep_trembl = info_cols[28]
			vep_uniparc = info_cols[29]
			vep_sift = info_cols[30]
			vep_domains = info_cols[31]
			vep_hgvs_offset = info_cols[32]
			# Create counter for each germline and tumor bam file with variant at sufficient VAF and set to zero
			g_sample_var_count = 0
			t_sample_var_count = 0
			# Create coverage variables
			same_tumor_cov = "no"
			gleich_germline_cov = "no"
			# Reset all variables
			g_bird = None
			germline_bam = None
			g_mpu_bases = None
			g_mpu_depth = None
			g_VAC = None
			g_VAF = None
			g_samtools_cmd = None
			g_samtools_proc = None
			g_out = None
			g_err = None
			g_mpu_out = None
			g_mpu = None
			same_bird = None
			same_mpu_chr = None
			same_mpu_pos = None
			same_mpu_ref = None
			same_mpu_depth = None
			same_mpu_bases = None
			same_VAC = None
			same_VAF = None
			g_mpu_chr = None
			g_mpu_pos = None
			g_mpu_ref = None
			g_mpu_depth = None
			g_mpu_bases = None
			g_VAC = None
			g_VAF = None
			t_bird = None
			tumor_bam = None
			t_samtools_cmd = None
			t_samtools_proc = None
			t_out = None
			t_err = None
			t_mpu_out = None
			t_mpu = None
			gleich_bird = None
			gleich_mpu_chr = None
			gleich_mpu_pos = None
			gleich_mpu_ref = None
			gleich_mpu_depth = None
			gleich_mpu_bases = None
			gleich_VAC = None
			gleich_VAF = None
			t_mpu_chr = None
			t_mpu_pos = None
			t_mpu_ref = None
			t_mpu_depth = None
			t_mpu_bases = None
			t_VAC = None
			t_VAF = None
			tumor_samples = None
			tumor_VAC = None
			tumor_VAF = None
			gleich_tumor_status = 'no'
			same_germline_status = 'no'
			tumor_in_germline_out = 'no'
			if vep_impact == 'MODERATE' or vep_impact == 'HIGH':
				# Search each tumor and germline bam file for the variant
				for g_bird, germline_bam in germline_bam_file.items():
					# Set all counting variables to zero
					g_mpu_bases = ''
					g_mpu_depth = 0
					g_VAC = 0
					g_VAF = 0
					# Put the command in a variable
					g_samtools_cmd = 'samtools mpileup -r ' + vep_chr+':'+vep_pos+'-'+vep_pos+' '+germline_bam
					# Use subprocess.Popen to ellicit shell commands 
					g_samtools_proc = subprocess.Popen([g_samtools_cmd], stdout=subprocess.PIPE, shell=True)
					# Use communicate to capture the output in a 'bytes' object
					(g_out, g_err) = g_samtools_proc.communicate()
					# Decode the 'bytes' object to a string
					g_mpu_out = g_out.decode("utf-8")
					g_mpu = g_mpu_out.rstrip()
					if g_bird[0:3] == vep_sample[0:3] and g_mpu != '':
						#print(g_mpu)
						same_bird = vep_sample
						same_mpu_chr = g_mpu.split('\t')[0]
						same_mpu_pos = g_mpu.split('\t')[1]
						same_mpu_ref = g_mpu.split('\t')[2]
						same_mpu_depth = int(g_mpu.split('\t')[3])
						if same_mpu_depth == 0:
							pass
						else:	
							same_mpu_bases = g_mpu.split('\t')[4].upper()
							same_VAC = same_mpu_bases.count(vep_alt)
							same_VAF = same_VAC/same_mpu_depth
							if same_VAF >= 0.05 and same_mpu_depth >= 4:
								g_sample_var_count = g_sample_var_count + 1
								same_germline_status = 'yes'
								#print('SAME GERMLINE STATUS WORKS: ' + same_germline_status + '\n' + '\n')
							if same_mpu_depth >= 4:
								same_tumor_cov = "yes"
							else:
								same_tumor_cov = "no"
					elif g_mpu == '' or int(g_mpu.split('\t')[3]) == 0:
						pass
					else:
						#print(g_mpu)
						g_mpu_chr = g_mpu.split('\t')[0]
						g_mpu_pos = g_mpu.split('\t')[1]
						g_mpu_ref = g_mpu.split('\t')[2]
						g_mpu_depth = int(g_mpu.split('\t')[3])
						g_mpu_bases = g_mpu.split('\t')[4].upper()
						g_VAC = g_mpu_bases.count(vep_alt)
						g_VAF = g_VAC/g_mpu_depth
						# Add to counter for each germline file with variant at sufficient VAF
						if g_VAF >= 0.10:
							g_sample_var_count = g_sample_var_count + 1
				# Search input tumor bam for each somatic called variant
				for t_bird, tumor_bam in tumor_bam_file.items():
					t_samtools_cmd = 'samtools mpileup -r ' + vep_chr+':'+vep_pos+'-'+vep_pos+' '+tumor_bam
					t_samtools_proc = subprocess.Popen([t_samtools_cmd], stdout=subprocess.PIPE, shell=True)
					(t_out, t_err) = t_samtools_proc.communicate()
					t_mpu_out = t_out.decode("utf-8")
					t_mpu = t_mpu_out.rstrip()
					if t_bird == vep_sample and t_mpu != '':
						gleich_bird = vep_sample
						gleich_mpu_chr = t_mpu.split('\t')[0]
						gleich_mpu_pos = t_mpu.split('\t')[1]
						gleich_mpu_ref = t_mpu.split('\t')[2]
						gleich_mpu_depth = int(t_mpu.split('\t')[3])
						gleich_mpu_bases = t_mpu.split('\t')[4].upper()
						gleich_VAC = gleich_mpu_bases.count(vep_alt)
						gleich_VAF = gleich_VAC/gleich_mpu_depth
						if gleich_VAF >= 0.05 and gleich_mpu_depth >= 4:
							t_sample_var_count = t_sample_var_count + 1
							gleich_tumor_status = 'yes'
							if tumor_samples != None:
								tumor_samples = tumor_samples + '|' + t_bird
							elif tumor_samples == None:
								tumor_samples = t_bird
							if tumor_VAC != None:
								tumor_VAC = tumor_VAC + '|' + str(gleich_VAC)
							elif tumor_VAC == None:
								tumor_VAC = str(gleich_VAC)
							if tumor_VAF != None:
								tumor_VAF = tumor_VAF + '|' + str(gleich_VAF)[0:5]
							elif tumor_VAF == None:
								tumor_VAF = str(gleich_VAF)[0:5]
							#print('Gleich Tumor STATUS WORKS: ' + gleich_tumor_status + '\n' + '\n')
						if gleich_mpu_depth >= 4:
							gleich_germline_cov = "yes"
						else:
							gleich_germline_cov = "no"
					elif t_mpu == '' or int(t_mpu.split('\t')[3]) == 0:
						pass
					#else:
						#t_mpu_chr = t_mpu.split('\t')[0]
						#print('vep_snv: ' + vep_snv)
						#print('Variant Call Sample: ' + vep_sample)
						#print('t_bird: ' + t_bird)
						#print('tumor_bam: ' + tumor_bam)
						#print('t_mpu: ' + t_mpu)
						#print('t_mpu_chr: ' + t_mpu_chr)
						#t_mpu_pos = t_mpu.split('\t')[1]
						#t_mpu_ref = t_mpu.split('\t')[2]
						#t_mpu_depth = int(t_mpu.split('\t')[3])
						#t_mpu_bases = t_mpu.split('\t')[4].upper()
						#t_VAC = t_mpu_bases.count(vep_alt)
						#t_VAF = t_VAC/t_mpu_depth
						#if t_VAF >= 0.05:
							#t_sample_var_count = t_sample_var_count + 1
							#if tumor_samples != None:
								#tumor_samples = tumor_samples + '|' + t_bird
							#elif tumor_samples == None:
								#tumor_samples = t_bird
							#if tumor_VAC != None:
								#tumor_VAC = tumor_VAC + '|' + str(t_VAC)
							#elif tumor_VAC == None:
								#tumor_VAC = str(t_VAC)
							#if tumor_VAF != None:
								#tumor_VAF = tumor_VAF + '|' + str(t_VAF)[0:5]
							#elif tumor_VAF == None:
								#tumor_VAF = str(t_VAF)[0:5]
				if gleich_tumor_status == 'yes' and same_germline_status == 'no':
					tumor_in_germline_out = 'yes'
				else:
					tumor_in_germline_out = 'no'
				if vep_symbol == '':
					vep_symbol = 'NA'
				# Perform final filters, if found at relevant freq in tumors and not found at slightly high freq in germline
				if t_sample_var_count > 0 and g_sample_var_count <= 0 and gleich_germline_cov == 'yes' and same_tumor_cov == 'yes' and tumor_in_germline_out == 'yes':
					outfile.write(vep_chr + '\t' + vep_pos + '\t' + vep_ref + '\t' + vep_alt + '\t' + vep_cons + '\t' + vep_impact + '\t' + vep_symbol + '\t' + vep_geneid + '\t' + str(t_sample_var_count) + '\t' + tumor_samples + '\t' + tumor_VAC + '\t' + tumor_VAF + '\n')

#logging.debug('This message should go to the log file')
#logging.info('So should this')
#logging.warning('And this, too')
