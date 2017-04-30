cd /home/proj/MDW_genomics/steepale/pathway_analysis

Tbird=$1

python ./scripts/custom_filter_somatic_snvs.py \
'./data/somaticseq_vcf/'$Tbird'_somaticseq_snv_vep.vcf' \
'./data/all_nonsyn_snvs_'$Tbird'.txt' \
$Tbird

# Grab the uniq lines (also corresponds to uniq snvs in this case) and sort
(grep "^#" './data/all_nonsyn_snvs_'$Tbird'.txt'; \
grep -v "^#" './data/all_nonsyn_snvs_'$Tbird'.txt' | sort | uniq) > \
'./data/all_nonsyn_snvs_'$Tbird'_final.txt'

# Remove redundant files
rm './data/all_nonsyn_snvs_'$Tbird'.txt'
