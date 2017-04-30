import sys
import os

#infile = open(sys.argv[1])
#outfile = open(sys.argv[2], 'w')

# Outfiles
#new_vep_file = {}
#open("./data/somaticseq_vcf/738-1_S1_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/741-1_S2_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/756-3_S3_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/766-1_S4_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/777-3_S14_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/787-2_S15_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/788-1_S16_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/794-1_S17_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/798-1_S5_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/833-1_S6_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/834-2_2_S12_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/834-2_S7_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/835-1_S18_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/841-3_S19_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/842-2_2_S25_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/842-2_S20_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/855-1_S8_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/863-1_S9_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/884-2_S21_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/901-2_2_S26_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/901-2_S22_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/906-1_S23_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/911-1_2_S13_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/911-1_S24_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/918-3_S10_somaticseq_snv_vep_labels.vcf", 'w')
#open("./data/somaticseq_vcf/927-2_S11_somaticseq_snv_vep_labels.vcf", 'w')


# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"

# Create a vep_vcf file dictionary
tumor_vep_file = {}
for tbird in open(tumor_birds_file):
	tbird = tbird.rstrip()
	tumor_vep_file[tbird] = './data/somaticseq_vcf/' + tbird + '_somaticseq_snv_vep.vcf'

# Create a dictionary of output files
new_vep_file = {}
for tbird in open(tumor_birds_file):
	tbird = tbird.rstrip()
	new_vep_file[tbird] = open('./data/somaticseq_vcf/'+tbird+'_somaticseq_snv_vep_labels.vcf', 'w')

for gbird in open(germline_birds_file):
	gbird = gbird.rstrip()
	for tbird in open(tumor_birds_file):
		tbird = tbird.rstrip()
		if gbird[0:3] == tbird[0:3]:
			#outfile = new_vep_file[tbird]
			#print(outfile)
			for vcf_vep_line in open(tumor_vep_file[tbird]):
				if vcf_vep_line.split('\t')[0] == '#CHROM':
					vep_vcf_line = vcf_vep_line.replace('NORMAL', gbird)
					(new_vep_file[tbird]).write(vep_vcf_line)
				else:
					(new_vep_file[tbird]).write(vcf_vep_line)
