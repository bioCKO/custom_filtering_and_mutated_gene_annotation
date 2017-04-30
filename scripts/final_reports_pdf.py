import sys
import re
import os
import subprocess
import numpy
# Modules for PDF production
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

# Input file
infile = sys.argv[1]

# Reference files:
###DONT_DELETE###chick_gene_file = "/home/users/a.steep/databases/ensembl/galgal5_ensembl_gene_transcript_protein.tsv"
SOMATIC_TEMP_chick_gene_file = "/home/users/a.steep/databases/ensembl/galgal5_ensembl_gene_transcript_protein_SOMATIC_TEMP.tsv"
###DONT_DELETE###ortholog_file = "/home/users/a.steep/databases/ensembl/chicken-human_orthologs.txt"
SOMATIC_TEMP_ortholog_file = "/home/users/a.steep/databases/ensembl/chicken-human_orthologs_SOMATIC_TEMP.txt"
###DONT_DELETE###ncbi2ensembl_file = "/home/users/a.steep/databases/ncbi/gene2ensembl"
SOMATIC_TEMP_ncbi2ensembl_file = "/home/users/a.steep/databases/ncbi/gene2ensembl_SOMATIC_TEMP.txt"
gene2refseq_file = "/home/users/a.steep/databases/ncbi/refseq/gene_RefSeqGene"
SOMATIC_TEMP_gene2refseq_file = open("/home/users/a.steep/databases/ncbi/refseq/gene_RefSeqGene_SOMATIC_TEMP.txt", 'w')
cosmic_cgc_file = "/home/users/a.steep/databases/cosmic/cosmic_CGC_gene_list.txt"
omin_file = "/home/users/a.steep/databases/omin/genemap2.txt"
sosnv_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final.txt"
soindel_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final.txt"
ss_vep_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final_vep.vcf"
si_vep_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final_vep.vcf"
birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
so_snv_and_indel_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_and_indels_final.txt"
somatic_lolliplots_dir = "/home/proj/MDW_genomics/steepale/GenVisR/analysis/plots/somatic_variants/"

refseq_files = [
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.1.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.2.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.3.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.4.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.5.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.6.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.7.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.8.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.9.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.10.genomic.gbff", \
"/home/users/a.steep/databases/ncbi/refseq/refseqgene.11.genomic.gbff"]

# Conditional variables
step1 = False
step2 = False

# Create empty dictionaries
ann_dict = {}
chick_gene_id2chick_gene_name = {}
chick_gene_id2chick_entrez_gene_id = {}
chick_gene_id2human_gene_id = {}
chick_gene_id2human_entrez_gene_id = {}
chick_gene_id2human_refseq_id = {}
chick_gene_id2human_refseq = {}
chick_gene_id2somatic_snv2vep = {}
chick_gene_id2somatic_indel2vep = {}

# Create a dictionary of vcf files
vcf_file = {}
vcf_file_indel = {}
# Iterate over samples labels (birds) and create dictionary of somatic snv vep files
# Iterate over samples labels (birds) and create dictionary of somatic snv vep files
for bird in open(birds_file):
    bird = bird.rstrip()
    #vcf_file.append(bird)
    vcf_file[bird] = "./data/somaticseq_vcf/" + bird + "_somaticseq_snv_vep.vcf"
    vcf_file_indel[bird] = "./data/somaticseq_vcf/" + bird + "_somaticseq_indel_vep.vcf"

# Empty variables
chick_gene_id = ''

# Output file path:
outfile_path = './results/mutated_gene_summaries/'

# Create a unique set of somatically mutated gene ids (ensembl) from infile
set_ens_gene_id = set()
# Iterate through lines of the infile
for in_line in open(infile):
    if in_line[0] != '#':
        in_cols = in_line.split('\t')
        in_ens_gene_id = in_cols[7]
        set_ens_gene_id.add(in_ens_gene_id)

# Build the appropriate dictionaries based on the input file
#for in_line in open(infile):
#    if in_line[0] != '#':
#        in_line = in_line.rstrip()
#        in_cols = in_line.split('\t')
#        in_gene_symbol = in_cols[6]
#        in_ensembl_gene_id = in_cols[7]
#        in_entrez_gene_id = in_cols[8]
#        print(len(ann_dict))
#        print(in_gene_symbol + '\t' + in_ensembl_gene_id)

# Iterate through the set of ensemble gene ids from input file
for chick_gene_id in set_ens_gene_id:
    print('dic1 ' + chick_gene_id)
    # Start building the annotation dictionaries for annotation purposes
    # Obtain chick central dogma of bio info
    for line in open(SOMATIC_TEMP_chick_gene_file):
        # Start on the correct line
        if line.split('\t')[0] != 'Gene_ID' and chick_gene_id == line.split('\t')[0]:
            line = line.rstrip()
            chick_gene_id = line.split('\t')[0]
            chick_transcript_id = line.split('\t')[1]
            chick_pro_id = line.split('\t')[2]
            chick_gene_name = line.split('\t')[3]
            chick_transcript_count = line.split('\t')[4]
            chick_gene_id2chick_gene_name[chick_gene_id] = chick_gene_name

# Save the generated dictionary
#numpy.save('./data/dictionaries/chick_gene_id2chick_gene_name.npy', chick_gene_id2chick_gene_name)
# Load the generated dictionary
#chick_gene_id2chick_gene_name = numpy.load('./data/dictionaries/chick_gene_id2chick_gene_name.npy').item()
#print(chick_gene_id2chick_gene_name)

# Iterate through the set of ensemble gene ids from input file
for chick_gene_id in set_ens_gene_id:
    print('dic2 ' + chick_gene_id)
    # Obtain ensembl ortholog gene name and ID
    for orth_line in open(SOMATIC_TEMP_ortholog_file):
        if orth_line.split('\t')[1] == chick_gene_id and orth_line.split('\t')[4].split('\n')[0] == "ortholog_one2one":
            human_gene_id = orth_line.split('\t')[3]
            human_gene_name = orth_line.split('\t')[2]
            # Create an empty dictionary with key value as ensembl gene id
            chick_gene_id2human_gene_id[chick_gene_id] = human_gene_id

for chick_gene_id in set_ens_gene_id:
    print('dic3 ' + chick_gene_id)
    # Obtain Entrez ID
    for ncbi2ensembl_line in open(SOMATIC_TEMP_ncbi2ensembl_file):
        ncbi2ensembl_line = ncbi2ensembl_line.rstrip()
        if ncbi2ensembl_line.split('\t')[2] == chick_gene_id:
            chick_entrez_gene_id = ncbi2ensembl_line.split('\t')[1]
            chick_gene_id2chick_entrez_gene_id[chick_gene_id] = chick_entrez_gene_id
            #print(chick_gene_id2chick_entrez_gene_id)
        if chick_gene_id in chick_gene_id2human_gene_id.keys():
            if ncbi2ensembl_line.split('\t')[2] == str(chick_gene_id2human_gene_id[chick_gene_id]):
                human_entrez_gene_id = ncbi2ensembl_line.split('\t')[1]
                chick_gene_id2human_entrez_gene_id[chick_gene_id] = human_entrez_gene_id
                #print(chick_gene_id2chick_entrez_gene_id)

# Load the generated dictionary
chick_gene_id2human_refseq = numpy.load('./data/dictionaries/chick_gene_id2human_refseq_2.npy').item()

# Iterate through the set of ensemble gene ids from input file
#for chick_gene_id in set_ens_gene_id:
#    print('refseq ' + chick_gene_id)
#    if chick_gene_id in chick_gene_id2chick_gene_name.keys():
#        print(chick_gene_id2chick_gene_name[chick_gene_id])
#    # Obtain human refseq info
#    for gene2refseq_line in open(gene2refseq_file):
#        gene2refseq_line = gene2refseq_line.rstrip()
#        # Good spot for empty variables so script doesn't hang
#        refseq1 = ''
#        refseq2 = ''
#        refseq3 = ''
#        refseq4 = ''
#        human_refseq = ''
#        if chick_gene_id in chick_gene_id2human_entrez_gene_id.keys():
#            if gene2refseq_line.split('\t')[1] == chick_gene_id2human_entrez_gene_id[chick_gene_id]:
#                human_refseq_id = gene2refseq_line.split('\t')[3]
#                chick_gene_id2human_refseq_id[chick_gene_id] = human_refseq_id
#                # Start building output value for dictionary
#                # Obtain a summary of the gene in refseq
#                ref_seq_cmd = "zgrep -w " + chick_gene_id2human_refseq_id[chick_gene_id] + " " + "/home/users/a.steep/databases/ncbi/refseq/refseqgene.*.genomic.gbff.gz" + " | cut -d':' -f1 | sort | uniq | sed 's/.gz//'"
#                #refseq_file = os.system(ref_seq_cmd)
#                proc = subprocess.Popen([ref_seq_cmd], stdout=subprocess.PIPE, shell=True) 
#                (out, err) = proc.communicate()
#                out = str(out)
#                refseq_file = out.replace("b'", "")[:-3]
#                print(refseq_file)
#                with open(refseq_file) as regseq_file_mem:
#                    while True:
#                        refseq_lines = regseq_file_mem.readlines(22000000)
#                        if not refseq_lines:
#                            break
#                        for refseq_line in refseq_lines:
#                            refseq_line = refseq_line.rstrip()
#                            if re.search("^LOCUS", refseq_line) and re.search(human_refseq_id.split('.')[0], refseq_line):
#                                step1 = True
#                            elif re.search("^PRIMARY", refseq_line):
#                                step1 = False
#                            elif step1:
#                                if re.search("DEFINITION", refseq_line):
#                                    refseq1 = refseq_line + '\n'
#                                elif re.search("SOURCE", refseq_line):
#                                    refseq2 = refseq_line + '\n' + '\n'
#                                elif re.search("Summary:", refseq_line):
#                                    refseq3 = refseq_line
#                                    step2 = True
#                                elif re.search("PRIMARY", refseq_line):
#                                    step2 = False
#                                elif step2:
#                                    refseq4 = refseq4 + refseq_line
#                                else:
#                                    step2 = False
#                human_refseq = refseq1 + refseq2 + refseq3 + refseq4.replace('   ', '')
#                if human_refseq != '':
#                    #ann_dict[chick_gene_id] = {chick_gene_name: {human_gene_id: {human_entrez_gene_id: {human_refseq_id: human_refseq}}}}
#                    chick_gene_id2human_refseq[chick_gene_id] = human_refseq
#                    print(chick_gene_id2human_refseq[chick_gene_id])
#                # Save the generated dictionary
#                #numpy.save('./data/dictionaries/chick_gene_id2human_refseq_3.npy', chick_gene_id2human_refseq)
#                # Load the generated dictionary
#               chick_gene_id2human_refseq = numpy.load('./data/dictionaries/chick_gene_id2human_refseq_2.npy').item()
                

# Start building the output file
# Story line for PDF output
Story=[]

# Loop though dictionary values and grab all relevant chick ensembl gene ids
# Get other chick and othrologue human gene identifiers
for chick_gene_id in set_ens_gene_id:
    if chick_gene_id in chick_gene_id2chick_gene_name.keys():
        chick_gene_name = str(chick_gene_id2chick_gene_name[chick_gene_id])
    else:
        chick_gene_name = ''
    print(chick_gene_id)
    print(chick_gene_name)
    if chick_gene_id in chick_gene_id2chick_entrez_gene_id.keys():
        chick_entrez_gene_id = str(chick_gene_id2chick_entrez_gene_id[chick_gene_id])
    else:
        chick_entrez_gene_id = ''
    if chick_gene_id in chick_gene_id2human_gene_id.keys():
        human_ensembl_gene_id = str(chick_gene_id2human_gene_id[chick_gene_id])
    else:
        human_ensembl_gene_id = ''
    if chick_gene_id in chick_gene_id2human_entrez_gene_id.keys():
        human_entrez_gene_id = str(chick_gene_id2human_entrez_gene_id[chick_gene_id])
    else:
        human_entrez_gene_id = ''
    # If there is no chick gene name associated with chick ensemble ID, use ensemble ID in output file

    if len(chick_gene_id2chick_gene_name[chick_gene_id]) > 0:
        #outfile = open(outfile_path + str(chick_gene_id2chick_gene_name[chick_gene_id]) + '_annotations.txt', 'w')
        outfile = SimpleDocTemplate(outfile_path + str(chick_gene_id2chick_gene_name[chick_gene_id]) + "_annotations.pdf",pagesize=letter,
                        rightMargin=40,leftMargin=40,
                        topMargin=40,bottomMargin=18)
    # If there is a chick gene symbol associated with chicken gene symbol, use gene symbol in output file name
    elif len(chick_gene_id2chick_gene_name[chick_gene_id]) == 0:
        #outfile = open(outfile_path + chick_gene_id + '_annotations.txt', 'w')
        outfile = SimpleDocTemplate(outfile_path + chick_gene_id + "_annotations.pdf",pagesize=letter,
                        rightMargin=40,leftMargin=40,
                        topMargin=40,bottomMargin=18)
    # Generate a style format (from imported module) for the output PDF
    styles=getSampleStyleSheet()
    # Further adjust the style guide
    styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

    # Capture Ensembl chick gene ID for output
    out_chick_ensembl_gene_id = 'Ensembl Chicken Gene ID: ' + chick_gene_id
    # Create the font size and add appropriate value in output
    pdf_text = '<font size=8>%s</font>' % out_chick_ensembl_gene_id 
    # Append value to storyline of PDF to be created
    Story.append(Paragraph(pdf_text, styles["Normal"]))

    # Capture Ensembl chick gene symbol for output
    if len(chick_gene_id2chick_gene_name[chick_gene_id]) > 0:
        out_ens_chick_gene_symbol = 'Ensembl Chicken Gene Name: ' + chick_gene_id2chick_gene_name[chick_gene_id]

        pdf_text = '<font size=8>%s</font>' % out_ens_chick_gene_symbol
        Story.append(Paragraph(pdf_text, styles["Normal"]))
        # Create a spacer line in PDF
        Story.append(Spacer(1, 12))

    # Capture human refseq gene summary for output
    if chick_gene_id in chick_gene_id2human_refseq.keys():
        # RefSeq title
        pdf_text = '<font size=8>RefSeq Human Ortholog Summary:</font>'
        Story.append(Paragraph(pdf_text, styles["Normal"]))

        out_refseq_human_gene = chick_gene_id2human_refseq[chick_gene_id]
        pdf_text = '<font size=8>%s</font>' % out_refseq_human_gene
        Story.append(Paragraph(pdf_text, styles["Normal"]))
        Story.append(Spacer(1, 12))
    
    # Grab COSMIC CGC annotation:
    # for chick_gene_id, values in ann_dict.items():
    #chick_gene_id = str(chick_gene_id) 
    #human_entrez_gene_id = str(chick_gene_id2human_entrez_gene_id[chick_gene_id])
    # Loop through COMSIC file
    for cgc_line in open(cosmic_cgc_file):
        cgc_line = cgc_line.rstrip()
        # If human entrez gene id from orthologue matches, then capture CGC info
        if cgc_line.split('\t')[2] == human_entrez_gene_id:
            cgc_cols = cgc_line.split('\t')
            cgc_symbol = cgc_cols[0]
            cgc_synonyms = cgc_cols[17]
            cgc_gene_name = cgc_cols[1]
            cgc_somatic = cgc_cols[5]
            cgc_germline = cgc_cols[6]
            cgc_som_types = cgc_cols[7]
            cgc_germ_types = cgc_cols[8]
            cgc_cancer_syndrome = cgc_cols[9]
            cgc_tissue = cgc_cols[10]
            cgc_genetics = cgc_cols[11]
            cgc_role = cgc_cols[12]
            cgc_mut_types = cgc_cols[13]
            cgc_trans_part = cgc_cols[14]
            cgc_other_germ = cgc_cols[15]
            cgc_other_syn = cgc_cols[16]
            # Output the CGC info:
            # CGC Title
            pdf_text = '<font size=8>COSMIC Cancer Gene Consensus:</font>' 
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Gene Symbol
            pdf_text = '<font size=8>%s</font>' % cgc_symbol
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Gene Name
            pdf_text = '<font size=8>%s</font>' % cgc_gene_name
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Gene Symbol Synonyms
            out_cgc_synonyms = 'Synonyms: ' + cgc_synonyms
            pdf_text = '<font size=8>%s</font>' % out_cgc_synonyms
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Cancer Types
            if cgc_somatic == 'yes':
                out_cgc_som_types = 'Cancers Associated with Somatic Mutations: ' + cgc_som_types
                pdf_text = '<font size=8>%s</font>' % out_cgc_som_types
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Germline Mutation Types
            if cgc_germline == 'yes':
                out_cgc_germ_types = 'GCancers Associated with Germline Mutations: ' + cgc_germ_types
                pdf_text = '<font size=8>%s</font>' % out_cgc_germ_types
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Known Cancer Syndromes
            if cgc_cancer_syndrome != '':
                out_cgc_germline = 'Known Cancer Syndromes: ' + cgc_germline
                pdf_text = '<font size=8>%s</font>' % out_cgc_germline
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Tissue Types
            if cgc_tissue != '':
                out_cgc_tissue = 'Tissue Types: ' + cgc_tissue
                pdf_text = '<font size=8>%s</font>' % out_cgc_tissue
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Role of Mutated Gene in Cancer
            if cgc_role != '':
                out_cgc_role = 'Role of Mutated Gene in Cancer: ' + cgc_role
                pdf_text = '<font size=8>%s</font>' % out_cgc_role
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Molecular Genetics
            if cgc_genetics != '':
                out_cgc_genetics = 'Molecular Genetics: ' + cgc_genetics
                pdf_text = '<font size=8>%s</font>' % out_cgc_genetics
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Mutation Types
            if cgc_mut_types != '':
                out_cgc_mut_types = 'Mutation Types: ' + cgc_mut_types
                pdf_text = '<font size=8>%s</font>' % out_cgc_mut_types
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Other Germline Mutations
            if cgc_other_germ == 'yes':
                out_cgc_other_germ = 'Other Germline Mutations: ' + cgc_other_germ
                pdf_text = '<font size=8>%s</font>' % out_cgc_other_germ
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Other Syndrome
            if cgc_other_syn != '':
                out_cgc_other_syn = 'Other Syndrome: ' + cgc_other_syn
                pdf_text = '<font size=8>%s</font>' % out_cgc_other_syn
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # Spacer
            Story.append(Spacer(1, 12))
    
    # Grab OMIN annotation:
    # for chick_gene_id, values in ann_dict.items():
    #chick_gene_id = str(chick_gene_id)
    # Lopp through OMIM File
    for omin_line in open(omin_file):
        if omin_line[0] != '#':
            omin_line = omin_line.rstrip()
            # If the orthologus human entrez or ensembl gene ID matches gene, collect OMIM variables
            if human_ensembl_gene_id != '' and re.search(human_ensembl_gene_id, omin_line):
                omin_cols = omin_line.split('\t')
                try:
                    omin_symbol = omin_cols[8]
                except IndexError:
                    omin_symbol = 'NA'
                try:
                    omin_com = omin_cols[11]
                except IndexError:
                    omin_com = 'NA'
                try:
                    omin_pheno = omin_cols[12]
                except IndexError:
                    omin_pheno = 'NA'
                #Output the OMIM info:
                # OMIM Title
                pdf_text = '<font size=8>Online Mendelian Inheritance in Man (OMIM):</font>' 
                Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM Gene Symbol
                pdf_text = '<font size=8>%s</font>' % omin_symbol
                Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM General Comments
                if omin_com != '':
                    out_omin_com = 'General Comments: ' + omin_com
                    pdf_text = '<font size=8>%s</font>' % out_omin_com
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM Mutation Associated Phenotypes
                if omin_pheno != '':
                    out_omin_pheno = 'Phenotypes Associated with Mutated Gene: ' + omin_pheno
                    pdf_text = '<font size=8>%s</font>' % out_omin_pheno
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                # Spacer
                Story.append(Spacer(1, 12))
            elif human_entrez_gene_id != '' and re.search('\t'+human_entrez_gene_id+'\t', omin_line):
                omin_cols = omin_line.split('\t')
                try:
                    omin_symbol = omin_cols[8]
                except IndexError:
                    omin_symbol = 'NA'
                try:
                    omin_com = omin_cols[11]
                except IndexError:
                    omin_com = 'NA'
                try:
                    omin_pheno = omin_cols[12]
                except IndexError:
                    omin_pheno = 'NA'
                #Output the OMIM info:
                # OMIM Title
                pdf_text = '<font size=8>Online Mendelian Inheritance in Man (OMIM):</font>' 
                Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM Gene Symbol
                pdf_text = '<font size=8>%s</font>' % omin_symbol
                Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM General Comments
                if omin_com != '':
                    out_omin_com = 'General Comments: ' + omin_com
                    pdf_text = '<font size=8>%s</font>' % out_omin_com
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                # OMIM Mutation Associated Phenotypes
                if omin_pheno != '':
                    out_omin_pheno = 'Phenotypes Associated with Mutated Gene: ' + omin_pheno
                    pdf_text = '<font size=8>%s</font>' % out_omin_pheno
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                # Spacer
                Story.append(Spacer(1, 12))
    
    # Annotate the Lolliplots of Somaticly Mutated Genes
    # Tutorial for reportlab: http://www.blog.pythonlibrary.org/2010/03/08/a-simple-step-by-step-reportlab-tutorial/
    # Documentation for reportlab: https://www.reportlab.com/docs/reportlab-userguide.pdf
    # Lilliplot Title
    if chick_gene_name == '':
        #print('GENE_NAME = 0')
        if re.search(chick_gene_id, str(os.listdir(somatic_lolliplots_dir))):
            pdf_text = '<font size=10>Lolliplot of Somatic Mutations in Gene Protein Products:</font>' 
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            Story.append(Spacer(1, 12))
    elif chick_gene_name != '':
        #print('GENE_NAME != 0')
        if re.search(chick_gene_name, str(os.listdir(somatic_lolliplots_dir))):
            pdf_text = '<font size=10>Lolliplot of Somatic Mutations in Gene Protein Products:</font>' 
            Story.append(Paragraph(pdf_text, styles["Normal"]))
            Story.append(Spacer(1, 12))

    # Iterate over lolliplot files
    for lolliplot in os.listdir(somatic_lolliplots_dir):
        # Query for files of appropriate gene lolliplot based on ensembl gene id or gene symbol
        if chick_gene_name == '':
            if lolliplot.startswith(chick_gene_id):
                #print('Inquiry of gene id')
                #print('Gene_id: ' + chick_gene_id)
                # Output the lolliplot
                out_lolli = Image(str(somatic_lolliplots_dir + lolliplot), width=5*inch, height=3*inch)
                Story.append(out_lolli)
        elif chick_gene_name != '':
            if lolliplot.startswith(chick_gene_name):
                #print('Inquiry of gene name')
                #print('Gene_id: ' + chick_gene_id)
                #print('Gene name: ' + chick_gene_name)
                # Output the lolliplot
                out_lolli = Image(str(somatic_lolliplots_dir + lolliplot), width=5*inch, height=3*inch)
                Story.append(out_lolli)

    # Output the somatic SNVs and INDELs for reference
    if re.search(chick_gene_name, str(os.listdir(somatic_lolliplots_dir))) or re.search(chick_gene_id, str(os.listdir(somatic_lolliplots_dir))):
        # Reset a somatic variant output variable before looping though somatic variants file
        out_somatic_vars = []
        n = 0
        # Loop through somatic SNV and INDEL final output file
        for so_snv_and_indel_line in open(so_snv_and_indel_file):
            if so_snv_and_indel_line[0] == '#':
                so_snv_and_indel_header = so_snv_and_indel_line
            elif so_snv_and_indel_line[0] != '#':
                so_snv_and_indel_line = so_snv_and_indel_line.rstrip()
                so_snv_and_indel_cols = so_snv_and_indel_line.split('\t')
                so_snv_and_indel_chick_gene_symbol = so_snv_and_indel_cols[6]
                so_snv_and_indel_chick_ensembl_gene_id = so_snv_and_indel_cols[7]
                # If the ensembl gene id matches in file, print variant line
                if so_snv_and_indel_chick_ensembl_gene_id == chick_gene_id:
                    n = n + 1
                    m = n - 1
                    out_somatic_vars.append(m)
                    out_somatic_vars[m] = so_snv_and_indel_line
        
        # Print the Header
        pdf_text = '<font size=10>%s</font>' % so_snv_and_indel_header
        Story.append(Paragraph(pdf_text, styles["Normal"]))
        # Print the Individual Variants
        for v in range(n):
            pdf_text = '<font size=6>%s</font>' % out_somatic_vars[v]
            Story.append(Paragraph(pdf_text, styles["Normal"]))
        Story.append(Spacer(1, 12))

    # Annotate the Somatic SNV Calls followed by VEP annotations
    # Print the Title
    pdf_text = '<font size=10>Detailed Annotation of Somatic SNVs:</font>' 
    Story.append(Paragraph(pdf_text, styles["Normal"]))
    Story.append(Spacer(1, 12))

    # Loop through the somatic SNV file
    for sosnv_line in open(sosnv_file):
        if sosnv_line[0] == '#':
            so_snv_line_header = sosnv_line
        if sosnv_line[0] != '#':
            sosnv_line = sosnv_line.rstrip()
            sosnv_cols = sosnv_line.split('\t')
            sosnv_chick_gene_symbol = sosnv_cols[6]
            sosnv_chick_ensembl_gene_id = sosnv_cols[7]
            # If the ensembl chick gene ID's match, collect info
            if sosnv_chick_ensembl_gene_id == chick_gene_id:
                sosnv_chr = sosnv_cols[0]
                sosnv_pos = sosnv_cols[1]
                sosnv_ref = sosnv_cols[2]
                sosnv_alt = sosnv_cols[3]
                # Capture unique snv info in variable
                sosnv_snv = sosnv_chr + sosnv_pos + sosnv_ref + sosnv_alt
                # Capture unique sample info in variable
                sosnv_sample = sosnv_cols[9]
                sosnv_sample_num = sosnv_sample.count('|') + 1
                # Note: Range treats numbers as zero-based
                sample = []
                # Collect info on each sample of variant
                for n in range(sosnv_sample_num):
                    sample.append(n)
                    sample[n] = sosnv_sample.split('|')[n]
                    # Add additional annotation from vep files
                    ss_vep_file2read = vcf_file[sample[n]]
                    # Loop through the VEP annotated file that is associated with captured sample above
                    for ss_vep_line in open(ss_vep_file2read):
                        if ss_vep_line[0] != '#': 
                            ss_vep_cols = ss_vep_line.split('\t')
                            ss_vep_chr = ss_vep_cols[0]
                            ss_vep_pos = ss_vep_cols[1]
                            ss_vep_ref = ss_vep_cols[3]
                            ss_vep_alt = ss_vep_cols[4]
                            ss_vep_snv = ss_vep_chr + ss_vep_pos + ss_vep_ref + ss_vep_alt
                            # If the variant from the somatic SNV file matches the variant from the sample's vep file,
                            # then collect vep annotation
                            if sosnv_snv == ss_vep_snv:
                                ss_vep_info = ss_vep_cols[7]
                                ss_vep_format = ss_vep_cols[8]
                                ss_vep_normal = ss_vep_cols[9]
                                ss_vep_tumor = ss_vep_cols[10]
                                if re.search('SOMATIC', ss_vep_info.split(';')[0]):
                                    ss_vep_somat = ss_vep_info.split(';')[0]
                                else:
                                    ss_vep_somat = 'NA'
                                if re.search('MVJSDU', ss_vep_info.split(';')[0]):
                                    ss_vep_tools = ss_vep_info.split(';')[0]
                                elif re.search('MVJSDU', ss_vep_info.split(';')[1]):
                                    ss_vep_tools = ss_vep_info.split(';')[1]
                                if re.search('NUM_TOOLS', ss_vep_info.split(';')[1]):
                                    ss_vep_tool_num = ss_vep_info.split(';')[1]
                                elif re.search('NUM_TOOLS', ss_vep_info.split(';')[2]):
                                    ss_vep_tool_num = ss_vep_info.split(';')[2]
                                if re.search('CSQ=', ss_vep_info.split(';')[2]):
                                    ss_vep_info_ann = ss_vep_info.split(';')[2].split('SQ=')[1]
                                    ss_vep_info_ann_num = ss_vep_info_ann.count(',') + 1
                                    info_ann = []
                                elif re.search('CSQ=', ss_vep_info.split(';')[3]):
                                    ss_vep_info_ann = ss_vep_info.split(';')[3].split('SQ=')[1]
                                    ss_vep_info_ann_num = ss_vep_info_ann.count(',') + 1
                                    info_ann = []
                                for n in range(ss_vep_info_ann_num):
                                    info_ann.append(n)
                                    info_ann[n] = ss_vep_info_ann.split(',')[n]
                                    info_ann2read = info_ann[n]
                                    info_cols = info_ann2read.split('|')
                                    ss_vep_allele = info_cols[0]
                                    ss_vep_cons = info_cols[1]
                                    ss_vep_impact = info_cols[2]
                                    ss_vep_symbol = info_cols[3]
                                    ss_vep_geneid = info_cols[4]
                                    ss_vep_feat_type = info_cols[5]
                                    ss_vep_feature = info_cols[6]
                                    ss_vep_biotype = info_cols[7]
                                    ss_vep_exon = info_cols[8]
                                    ss_vep_intron = info_cols[9]
                                    ss_vep_HGVSc = info_cols[10]
                                    ss_vep_HGVSp = info_cols[11]
                                    ss_vep_cDNA_pos = info_cols[12]
                                    ss_vep_CDS_pos = info_cols[13]
                                    ss_vep_protein_pos = info_cols[14]
                                    ss_vep_aminos = info_cols[15]
                                    ss_vep_codons = info_cols[16]
                                    ss_vep_existing_var = info_cols[17]
                                    ss_vep_distance = info_cols[18]
                                    ss_vep_strand = info_cols[19]
                                    ss_vep_flags = info_cols[20]
                                    ss_vep_symbol_source = info_cols[21]
                                    ss_vep_HGNC_ID = info_cols[22]
                                    ss_vep_tsl = info_cols[23]
                                    ss_vep_appris = info_cols[24]
                                    ss_vep_ccds = info_cols[25]
                                    ss_vep_ensp = info_cols[26]
                                    ss_vep_swissprot = info_cols[27]
                                    ss_vep_trembl = info_cols[28]
                                    ss_vep_uniparc = info_cols[29]
                                    ss_vep_sift = info_cols[30]
                                    ss_vep_domains = info_cols[31]
                                    ss_vep_hgvs_offset = info_cols[32]
                                    # Reset these just in case
                                    so_snv2vep1 = ''
                                    so_snv2vep2 = ''
                                    so_snv2vep3 = ''
                                    so_snv2vep4 = ''
                                    so_snv2vep5 = ''
                                    so_snv2vep6 = ''
                                    so_snv2vep7 = ''
                                    so_snv2vep8 = ''
                                    so_snv2vep9 = ''
                                    so_snv2vep10 = ''
                                    so_snv2vep11 = ''
                                    so_snv2vep12 = ''
                                    so_snv2vep13 = ''
                                    so_snv2vep14 = ''
                                    so_snv2vep15 = ''
                                    # Only print annotation on non-synonymous variables
                                    if ss_vep_impact == 'MODERATE' or ss_vep_impact == 'HIGH':
                                        so_snv2vep1 = so_snv_line_header
                                        so_snv2vep2 = sosnv_line
                                        so_snv2vep3 = 'VEP Annotation'
                                        if ss_vep_tool_num != '':
                                            so_snv2vep4 = ss_vep_tool_num + ' (Number of variant callers)'
                                        else:
                                            so_snv2vep4 = ''
                                        if ss_vep_tools != '':
                                            so_snv2vep5 = ss_vep_tools + ' (One letter abbreviations for each somatic-snv-caller)'
                                        else:
                                            so_snv2vep5 = ''
                                        if ss_vep_somat == 'SOMATIC':
                                            so_snv2vep6 = 'SomaticSeq predicted this variant to be somatic'
                                        else:
                                            so_snv2vep6 = ''
                                        if ss_vep_strand != '':
                                            so_snv2vep7 = 'DNA strand: ' + ss_vep_strand
                                        else:
                                            so_snv2vep7 = ''
                                        if ss_vep_sift != '':
                                            so_snv2vep8 = 'SIFT Loss-of-function prediciton: ' + ss_vep_sift
                                        else:
                                            so_snv2vep8 = ''
                                        if ss_vep_feat_type != '':
                                            so_snv2vep9 = 'Feature type: ' + ss_vep_feat_type
                                        else:
                                            so_snv2vep9 = ''
                                        if ss_vep_feature != '':
                                            so_snv2vep10 = 'Feature: ' + ss_vep_feature
                                        else:
                                            so_snv2vep10 = ''
                                        if ss_vep_biotype != '':
                                            so_snv2vep11 = 'Biotype: ' + ss_vep_biotype
                                        else:
                                            so_snv2vep11 = ''
                                        if ss_vep_exon != '':
                                            so_snv2vep12 = 'Exon: ' + ss_vep_exon
                                        else:
                                            so_snv2vep12 = ''
                                        if ss_vep_ensp != '' or ss_vep_swissprot != '' or ss_vep_trembl != '' or ss_vep_uniparc != '':
                                            so_snv2vep13 = 'Estimated protein identifiers:' 
                                            so_snv2vep13_1 = 'Ensembl: ' + ss_vep_ensp
                                            so_snv2vep13_2 = 'SwissProt: ' + ss_vep_swissprot
                                            so_snv2vep13_3 = 'TREMBL: ' + ss_vep_trembl
                                            so_snv2vep13_4 = 'UNIPARC: ' + ss_vep_uniparc
                                        else:
                                            so_snv2vep13 = ''
                                            so_snv2vep13_1 = ''
                                            so_snv2vep13_2 = ''
                                            so_snv2vep13_3 = ''
                                            so_snv2vep13_4 = ''
                                        if ss_vep_HGVSp != '':
                                            so_snv2vep14 = 'Protein alteration (HGVS): ' + ss_vep_HGVSp
                                        else:
                                            so_snv2vep14 = ''
                                        if ss_vep_domains != '':
                                            so_snv2vep15 = 'Mutation in protein domains: ' + ss_vep_domains
                                        else:
                                            so_snv2vep15 = ''
                                        # Consolidate output values
                                        so_snv2vep_fin = so_snv2vep1 + '<br/>' + so_snv2vep2 + '<br/>' + '<br/>' + so_snv2vep3 + '<br/>' + so_snv2vep4 + '<br/>' + so_snv2vep5 + '<br/>' + so_snv2vep6 + '<br/>' + so_snv2vep7 + '<br/>' + so_snv2vep8 + '<br/>' + so_snv2vep9 + '<br/>' + so_snv2vep10 + '<br/>' + so_snv2vep11 + '<br/>' + so_snv2vep12 + '<br/>' + so_snv2vep13 + '<br/>' + so_snv2vep14 + '<br/>' + so_snv2vep15
                                        # Associate consolidated output values with chicken ensembl gene ID in dictionary to eliminate redundency
                                        chick_gene_id2somatic_snv2vep[chick_gene_id] = {sosnv_snv: so_snv2vep_fin}
                
                if sosnv_snv in chick_gene_id2somatic_snv2vep[chick_gene_id].keys():
                    # Print the vep annotations for each somatic variant
                    pdf_text = '<font size=7>%s<br/> \
                    </font>' % (str(chick_gene_id2somatic_snv2vep[chick_gene_id][sosnv_snv]))
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                    Story.append(Spacer(1, 12))

    # Annotate the Somatic Indel Calls followed by VEP annotations
    # Print the Title
    pdf_text = '<font size=10>Detailed Annotation of Somatic INDELs:</font>' 
    Story.append(Paragraph(pdf_text, styles["Normal"]))
    Story.append(Spacer(1, 12))
    # Loop through the somatic INDEL files
    for soindel_line in open(soindel_file):
        if soindel_line[0] == '#':
            so_indel_line_header = soindel_line
        if soindel_line[0] != '#':
            soindel_line = soindel_line.rstrip()
            soindel_cols = soindel_line.split('\t')
            soindel_chick_gene_symbol = soindel_cols[6]
            soindel_chick_ensembl_gene_id = soindel_cols[7]
            # If the ensembl chick gene ID's match, collect info
            if soindel_chick_ensembl_gene_id == chick_gene_id:
                soindel_chr = soindel_cols[0]
                soindel_pos = soindel_cols[1]
                soindel_ref = soindel_cols[2]
                soindel_alt = soindel_cols[3]
                # Capture unique indel info in variable
                soindel_indel = soindel_chr + soindel_pos + soindel_ref + soindel_alt
                # Capture unique sample info in variable
                soindel_sample = soindel_cols[9]
                soindel_sample_num = soindel_sample.count('|') + 1
                # Note: Range treats numbers as zero-based
                sample_indel = []
                # Collect info on each sample of variant
                for n in range(soindel_sample_num):
                    sample_indel.append(n)
                    sample_indel[n] = soindel_sample.split('|')[n]
                    # Add additional annotation from vep files
                    si_vep_file2read = vcf_file_indel[sample_indel[n]]
                    # Loop through the VEP annotated file that is associated with captured sample above
                    for si_vep_line in open(si_vep_file2read):
                        if si_vep_line[0] != '#': 
                            si_vep_cols = si_vep_line.split('\t')
                            si_vep_chr = si_vep_cols[0]
                            si_vep_pos = si_vep_cols[1]
                            si_vep_ref = si_vep_cols[3]
                            si_vep_alt = si_vep_cols[4]
                            si_vep_indel = si_vep_chr + si_vep_pos + si_vep_ref + si_vep_alt
                            if soindel_indel == si_vep_indel:
                                #print('soindel_indel ' + soindel_indel)
                                #print('si_vep_indel ' + si_vep_indel)
                                si_vep_info = si_vep_cols[7]
                                si_vep_format = si_vep_cols[8]
                                si_vep_normal = si_vep_cols[9]
                                si_vep_tumor = si_vep_cols[10]
                                if re.search('SOMATIC', si_vep_info.split(';')[0]):
                                    si_vep_somat = si_vep_info.split(';')[0]
                                else:
                                    si_vep_somat = 'NA'
                                if re.search('MVDL', si_vep_info.split(';')[0]):
                                    si_vep_tools = si_vep_info.split(';')[0]
                                elif re.search('MVDL', si_vep_info.split(';')[1]):
                                    si_vep_tools = si_vep_info.split(';')[1]
                                if re.search('NUM_TOOLS', si_vep_info.split(';')[1]):
                                    si_vep_tool_num = si_vep_info.split(';')[1]
                                elif re.search('NUM_TOOLS', si_vep_info.split(';')[2]):
                                    si_vep_tool_num = si_vep_info.split(';')[2]
                                if re.search('CSQ=', si_vep_info.split(';')[2]):
                                    si_vep_info_ann = si_vep_info.split(';')[2].split('SQ=')[1]
                                    si_vep_info_ann_num = si_vep_info_ann.count(',') + 1
                                    info_ann_indel = []
                                elif re.search('CSQ=', si_vep_info.split(';')[3]):
                                    si_vep_info_ann = si_vep_info.split(';')[3].split('SQ=')[1]
                                    si_vep_info_ann_num = si_vep_info_ann.count(',') + 1
                                    info_ann_indel = []
                                for n in range(si_vep_info_ann_num):
                                    info_ann_indel.append(n)
                                    info_ann_indel[n] = si_vep_info_ann.split(',')[n]
                                    info_ann2read = info_ann_indel[n]
                                    info_cols = info_ann2read.split('|')
                                    si_vep_allele = info_cols[0]
                                    si_vep_cons = info_cols[1]
                                    si_vep_impact = info_cols[2]
                                    si_vep_symbol = info_cols[3]
                                    si_vep_geneid = info_cols[4]
                                    si_vep_feat_type = info_cols[5]
                                    si_vep_feature = info_cols[6]
                                    si_vep_biotype = info_cols[7]
                                    si_vep_exon = info_cols[8]
                                    si_vep_intron = info_cols[9]
                                    si_vep_HGVSc = info_cols[10]
                                    si_vep_HGVSp = info_cols[11]
                                    si_vep_cDNA_pos = info_cols[12]
                                    si_vep_CDS_pos = info_cols[13]
                                    si_vep_protein_pos = info_cols[14]
                                    si_vep_aminos = info_cols[15]
                                    si_vep_codons = info_cols[16]
                                    si_vep_existing_var = info_cols[17]
                                    si_vep_distance = info_cols[18]
                                    si_vep_strand = info_cols[19]
                                    si_vep_flags = info_cols[20]
                                    si_vep_symbol_source = info_cols[21]
                                    si_vep_HGNC_ID = info_cols[22]
                                    si_vep_tsl = info_cols[23]
                                    si_vep_appris = info_cols[24]
                                    si_vep_ccds = info_cols[25]
                                    si_vep_ensp = info_cols[26]
                                    si_vep_swissprot = info_cols[27]
                                    si_vep_trembl = info_cols[28]
                                    si_vep_uniparc = info_cols[29]
                                    si_vep_sift = info_cols[30]
                                    si_vep_domains = info_cols[31]
                                    si_vep_hgvs_offset = info_cols[32]
                                    # Reset these just in case
                                    si_indel2vep1 = ''
                                    si_indel2vep2 = ''
                                    si_indel2vep3 = ''
                                    si_indel2vep4 = ''
                                    si_indel2vep5 = ''
                                    si_indel2vep6 = ''
                                    si_indel2vep7 = ''
                                    si_indel2vep8 = ''
                                    si_indel2vep9 = ''
                                    si_indel2vep10 = ''
                                    si_indel2vep11 = ''
                                    si_indel2vep12 = ''
                                    si_indel2vep13 = ''
                                    si_indel2vep14 = ''
                                    si_indel2vep15 = ''
                                    if si_vep_impact == 'MODERATE' or si_vep_impact == 'HIGH':
                                        si_indel2vep1 = so_indel_line_header
                                        si_indel2vep2 = soindel_line
                                        si_indel2vep3 = 'VEP Annotation'
                                        if si_vep_tool_num != '':
                                            si_indel2vep4 = si_vep_tool_num + ' (Number of variant callers)'
                                        else:
                                            si_indel2vep4 = ''
                                        if si_vep_tools != '':
                                            si_indel2vep5 = si_vep_tools + ' (One letter abbreviations for each somatic-indel-caller)'
                                        else:
                                            si_indel2vep5 = ''
                                        if si_vep_somat == 'SOMATIC':
                                            si_indel2vep6 = 'SomaticSeq predicted this variant to be somatic'
                                        else:
                                            si_indel2vep6 = ''
                                        if si_vep_strand != '':
                                            si_indel2vep7 = 'DNA strand: ' + si_vep_strand
                                        else:
                                            si_indel2vep7 = ''
                                        if si_vep_sift != '':
                                            si_indel2vep8 = 'SIFT Loss-of-function prediciton: ' + si_vep_sift
                                        else:
                                            si_indel2vep8 = ''
                                        if si_vep_feat_type != '':
                                            si_indel2vep9 = 'Feature type: ' + si_vep_feat_type
                                        else:
                                            si_indel2vep9 = ''
                                        if si_vep_feature != '':
                                            si_indel2vep10 = 'Feature: ' + si_vep_feature
                                        else:
                                            si_indel2vep10 = ''
                                        if si_vep_biotype != '':
                                            si_indel2vep11 = 'Biotype: ' + si_vep_biotype
                                        else:
                                            si_indel2vep11 = ''
                                        if si_vep_exon != '':
                                            si_indel2vep12 = 'Exon: ' + si_vep_exon
                                        else:
                                            si_indel2vep12 = ''
                                        if si_vep_ensp != '' or si_vep_swissprot != '' or si_vep_trembl != '' or si_vep_uniparc != '':
                                            si_indel2vep13 = 'Estimated protein identifiers:'
                                            si_indel2vep13_1 = 'Ensembl: ' + si_vep_ensp
                                            si_indel2vep13_2 = 'SwissProt: ' + si_vep_swissprot
                                            si_indel2vep13_3 = 'TREMBL: ' + si_vep_trembl
                                            si_indel2vep13_4 = 'UNIPARC: ' + si_vep_uniparc
                                        else:
                                            si_indel2vep13 = ''
                                            si_indel2vep13_1 = ''
                                            si_indel2vep13_2 = ''
                                            si_indel2vep13_3 = ''
                                            si_indel2vep13_4 = ''
                                        if si_vep_HGVSp != '':
                                            si_indel2vep14 = 'Protein alteration (HGVS): ' + si_vep_HGVSp
                                        else:
                                            si_indel2vep14 = ''
                                        if si_vep_domains != '':
                                            si_indel2vep15 = 'Mutation in protein domains: ' + si_vep_domains
                                        else:
                                            si_indel2vep15 = ''
                                        si_indel2vep_fin = si_indel2vep1 + '<br/>' + si_indel2vep2 + '<br/>' + '<br/>' + si_indel2vep3 + '<br/>' + si_indel2vep4 + '<br/>' + si_indel2vep5 + '<br/>' + si_indel2vep6 + '<br/>' + si_indel2vep7 + '<br/>' + si_indel2vep8 + '<br/>' + si_indel2vep9 + '<br/>' + si_indel2vep10 + '<br/>' + si_indel2vep11 + '<br/>' + si_indel2vep12 + '<br/>' + si_indel2vep13 + '<br/>' + si_indel2vep14 + '<br/>' + si_indel2vep15 
                                        chick_gene_id2somatic_indel2vep[chick_gene_id] = {soindel_indel: si_indel2vep_fin}

                # Write the somatic Indel informatation to output file
                if soindel_indel in chick_gene_id2somatic_indel2vep[chick_gene_id]:
                    # Print the vep annotations for each somatic variant
                    pdf_text = '<font size=7>%s<br/> \
                    </font>' % (str(chick_gene_id2somatic_indel2vep[chick_gene_id][soindel_indel])) 
                    Story.append(Paragraph(pdf_text, styles["Normal"]))
                    Story.append(Spacer(1, 12))
                
    # Write the PDF output file
    outfile.build(Story)
