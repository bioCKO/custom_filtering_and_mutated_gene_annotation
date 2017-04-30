cd /home/proj/MDW_genomics/steepale/pathway_analysis

java -Xmx2g -jar /home/users/a.steep/Apps/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
-R /home/proj/MDW_genomics/steepale/galgal5/galgal5.fa \
-I /home/proj/MDW_genomics/xu/final_bam/${1}_Bwa_RG_dedupped_realigned.bam \
--variant ./data/somaticseq_vcf/${1}_somaticseq_snv.vcf.gz \
-L ./data/somaticseq_vcf/${1}_somaticseq_snv.vcf.gz \
-o ./data/somaticseq_vcf/${1}_somaticseq_snv_phased.vcf.gz \
--phaseQualityThresh 20.0

qstat -f ${PBS_JOBID}
