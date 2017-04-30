import sys
import re
import os
import subprocess
# Modules for PDF production
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch


# Input file
infile = sys.argv[1]

# Reference files:
#chick_gene_file = "/home/users/a.steep/databases/ensembl/galgal5_ensembl_gene_transcript_protein.tsv"
chick_gene_file = "/home/users/a.steep/databases/ensembl/galgal5_ensembl_gene_transcript_protein.tsv"
ortholog_file = "/home/users/a.steep/databases/ensembl/chicken-human_orthologs.txt"
#ncbi2ensembl_file = "/home/users/a.steep/databases/ncbi/gene2ensembl"
ncbi2ensembl_file = "/home/users/a.steep/databases/ncbi/gene2ensembl_uniq.txt"
gene2refseq_file = "/home/users/a.steep/databases/ncbi/refseq/gene_RefSeqGene"
cosmic_cgc_file = "/home/users/a.steep/databases/cosmic/cosmic_CGC_gene_list.txt"
omin_file = "/home/users/a.steep/databases/omin/genemap2.txt"
sosnv_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final.txt"
soindel_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final.txt"
ss_vep_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final_vep.vcf"
si_vep_file = "/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final_vep.vcf"
birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"


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

# Build the appropriate dictionaries based on the input file
for in_line in open(infile):
    if in_line[0] != '#':
        in_line = in_line.rstrip()
        in_cols = in_line.split('\t')
        in_gene_symbol = in_cols[6]
        in_ensembl_gene_id = in_cols[7]
        in_entrez_gene_id = in_cols[8]
        print(len(ann_dict))
        print(in_gene_symbol + '\t' + in_ensembl_gene_id)
        # Start building the annotation dictionaries for annotation purposes
        for line in open(chick_gene_file):
            #print(ann_dict)
            # Start on the correct line
            if line.split('\t')[0] != 'Gene_ID' and in_ensembl_gene_id == line.split('\t')[0]:
                #print(line)
                line = line.rstrip()
                chick_gene_id = line.split('\t')[0]
                chick_transcript_id = line.split('\t')[1]
                chick_pro_id = line.split('\t')[2]
                chick_gene_name = line.split('\t')[3]
                chick_transcript_count = line.split('\t')[4]
                chick_gene_id2chick_gene_name[chick_gene_id] = chick_gene_name
                # Obtain ensembl ortholog gene name and ID
                for orth_line in open(ortholog_file):
                    if orth_line.split('\t')[1] == chick_gene_id and orth_line.split('\t')[4].split('\n')[0] == "ortholog_one2one":
                        human_gene_id = orth_line.split('\t')[3]
                        human_gene_name = orth_line.split('\t')[2]
                        # Create an empty dictionary with key value as ensembl gene id
                        chick_gene_id2human_gene_id[chick_gene_id] = human_gene_id
                        # Obtain Entrez ID
                        for ncbi2ensembl_line in open(ncbi2ensembl_file):
                            ncbi2ensembl_line = ncbi2ensembl_line.rstrip()
                            if ncbi2ensembl_line.split('\t')[1] == chick_gene_id:
                                chick_entrez_gene_id = ncbi2ensembl_line.split('\t')[0]
                                chick_gene_id2chick_entrez_gene_id[chick_gene_id] = chick_entrez_gene_id
                            elif ncbi2ensembl_line.split('\t')[1] == human_gene_id:
                                human_entrez_gene_id = ncbi2ensembl_line.split('\t')[0]
                                chick_gene_id2human_entrez_gene_id[chick_gene_id] = human_entrez_gene_id
                                for gene2refseq_line in open(gene2refseq_file):
                                    gene2refseq_line = gene2refseq_line.rstrip()
                                    # Good spot for empty variables so script doesn't hang
                                    refseq1 = ''
                                    refseq2 = ''
                                    refseq3 = ''
                                    refseq4 = ''
                                    human_refseq = ''
                                    if gene2refseq_line.split('\t')[1] == chick_gene_id2human_entrez_gene_id[chick_gene_id]:
                                        human_refseq_id = gene2refseq_line.split('\t')[3]
                                        chick_gene_id2human_refseq_id[chick_gene_id] = human_refseq_id
                                        # Start building output file
                                        # Obtain a summary of the gene in refseq
                                        ref_seq_cmd = "zgrep -w " + chick_gene_id2human_refseq_id[chick_gene_id] + " " + "/home/users/a.steep/databases/ncbi/refseq/refseqgene.*.genomic.gbff.gz" + " | cut -d':' -f1 | sort | uniq | sed 's/.gz//'"
                                        #refseq_file = os.system(ref_seq_cmd)
                                        proc = subprocess.Popen([ref_seq_cmd], stdout=subprocess.PIPE, shell=True) 
                                        (out, err) = proc.communicate()
                                        out = str(out)
                                        refseq_file = out.replace("b'", "")[:-3]
                                        print(refseq_file)
                                        with open(refseq_file) as regseq_file_mem:
                                            while True:
                                                refseq_lines = regseq_file_mem.readlines(22000000)
                                                if not refseq_lines:
                                                    break
                                                for refseq_line in refseq_lines:
                                                    refseq_line = refseq_line.rstrip()
                                                    if re.search("^LOCUS", refseq_line) and re.search(human_refseq_id.split('.')[0], refseq_line):
                                                        step1 = True
                                                    elif re.search("^PRIMARY", refseq_line):
                                                        step1 = False
                                                    elif step1:
                                                        if re.search("DEFINITION", refseq_line):
                                                            refseq1 = refseq_line + '\n'
                                                        elif re.search("SOURCE", refseq_line):
                                                            refseq2 = refseq_line + '\n' + '\n'
                                                        elif re.search("Summary:", refseq_line):
                                                            refseq3 = refseq_line
                                                            step2 = True
                                                        elif re.search("PRIMARY", refseq_line):
                                                            step2 = False
                                                        elif step2:
                                                            refseq4 = refseq4 + refseq_line
                                                    else:
                                                        step2 = False
                                        human_refseq = refseq1 + refseq2 + refseq3 + refseq4.replace('   ', '')
                                        if human_refseq != '':
                                            ann_dict[chick_gene_id] = {chick_gene_name: {human_gene_id: {human_entrez_gene_id: {human_refseq_id: human_refseq}}}}
                                            chick_gene_id2human_refseq[chick_gene_id] = human_refseq
                                        else:
                                            break


# Start building the output file
# Story line for PDF output
Story=[]

# Loop though dictionary values and grab all relevant chick ensembl gene ids
# Get other chick and othrologue human gene identifiers
for chick_gene_id, values in ann_dict.items():
    chick_gene_id = str(chick_gene_id)
    chick_entrez_gene_id = str(chick_gene_id2chick_entrez_gene_id[chick_gene_id])
    human_ensembl_gene_id = str(chick_gene_id2human_gene_id[chick_gene_id])
    human_entrez_gene_id = str(chick_gene_id2human_entrez_gene_id[chick_gene_id])
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

    # Capture Ensembl chick gene symbol fot output
    if len(chick_gene_id2chick_gene_name[chick_gene_id]) > 0:
        out_ens_chick_gene_symbol = 'Ensembl Chicken Gene Name: ' + chick_gene_id2chick_gene_name[chick_gene_id]

        pdf_text = '<font size=8>%s</font>' % out_ens_chick_gene_symbol
        Story.append(Paragraph(pdf_text, styles["Normal"]))
        # Create a spacer line in PDF
        Story.append(Spacer(1, 12))

    # Capture human refseq gene summary fot output
    if len(chick_gene_id2human_refseq[chick_gene_id]) > 0:
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
            pdf_text = '<font size=8>COSMIC Cancer Gene Consensus</font>' 
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
            # CGC Somatic Mutation Types
            if cgc_somatic == 'yes':
                out_cgc_som_types = 'Somatic Mutation Types: ' + cgc_som_types
                pdf_text = '<font size=8>%s</font>' % out_cgc_som_types
                Story.append(Paragraph(pdf_text, styles["Normal"]))
            # CGC Germline Mutation Types
            if cgc_germline == 'yes':
                out_cgc_germ_types = 'Germline Mutation Types: ' + cgc_germ_types
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
    
    # Grab OMIN annotation:
    # for chick_gene_id, values in ann_dict.items():
    #chick_gene_id = str(chick_gene_id)
    # Lopp through OMIM File
    for omin_line in open(omin_file):
        if omin_line[0] != '#':
            omin_line = omin_line.rstrip()
            # If the orthologus human entrez or ensembl gene ID matches gene, collect OMIM variables
            if re.search(human_ensembl_gene_id, omin_line) or re.search('\t'+human_entrez_gene_id+'\t', omin_line):
                omin_cols = omin_line.split('\t')
                omin_symbol = omin_cols[8]
                omin_com = omin_cols[11]
                omin_pheno = omin_cols[12]
                #Output the OMIM info:
                # OMIM Title
                pdf_text = '<font size=8>OMIM Annotation</font>' 
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
    
    # Write the PDF output file
    outfile.build(Story)
