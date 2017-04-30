cd /home/proj/MDW_genomics/steepale/pathway_analysis

# Annotate with VEP
perl /home/users/a.steep/Apps/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i "./data/somaticseq_vcf/"$1"_somaticseq_snv.vcf" \
-o "./data/somaticseq_vcf/"$1"_somaticseq_snv_vep.vcf" \
--vcf \
--cache \
--species gallus_gallus \
--force_overwrite \
--plugin ProteinSeqs,"./data/proteinseqs/"$1"_somatic_snvs_reference.fa","./data/proteinseqs/"$1"_somatic_snvs_mutated.fa" \
--protein \
--hgvs \
--domains \
--ccds \
--uniprot \
--tsl \
--appris \
--sift b
