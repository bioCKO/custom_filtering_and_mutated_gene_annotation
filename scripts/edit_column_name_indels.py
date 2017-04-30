import sys
import os

#infile = open(sys.argv[1])
#outfile = open(sys.argv[2], 'w')

# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"

# Create a vep_vcf file dictionary
tumor_vep_file = {}
for tbird in open(tumor_birds_file):
	tbird = tbird.rstrip()
	tumor_vep_file[tbird] = './data/somaticseq_vcf/' + tbird + '_somaticseq_indel_vep.vcf'

# Create a dictionary of output files
new_vep_file = {}
for tbird in open(tumor_birds_file):
	tbird = tbird.rstrip()
	new_vep_file[tbird] = open('./data/somaticseq_vcf/'+tbird+'_somaticseq_indel_vep_labels.vcf', 'w')

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
