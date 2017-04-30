rm indels_germline_samples_pass.txt
rm indels_tumor_samples_pass.txt
rm ./data/filter/*.int
rm ./data/filter/*.test

echo "Counts" > ./indels_tumor_counts.txt
echo "Counts" > ./indels_germline_counts.txt
echo "FAIL" > ./data/somatic_indelss_fail_germline_sample_freq.txt
# Attempt 2016/12/15
# Print the header for variants in pass and fail files
# Header for filtered variants that passed (pass file 1)
echo "#Somatic indels that passed custom filtering" > ./data/somatic_indels_pass.txt
echo "#Filters include:" >> ./data/somatic_indels_pass.txt
echo "#Somatic indels cannot appear in more than 8 germline samples" >> ./data/somatic_indels_pass.txt
echo "#Variant allele must have at least 2 read bases" >> ./data/somatic_indels_pass.txt
echo "#Variant allele frequency must be at least 10%" >> ./data/somatic_indels_pass.txt
# Header for filtered variants that passed (final pass file)
echo "#Somatic indels that passed custom filtering" > ./data/somatic_indels_pass_final.txt
echo "#Filters include:" >> ./data/somatic_indels_pass_final.txt
echo "#Variant cannot appear in more than 8 germline samples" >> ./data/somatic_indels_pass_final.txt
echo "#Variant allele must have at least 2 read bases" >> ./data/somatic_indels_pass_final.txt
echo "#Variant allele frequency must be at least 10%" >> ./data/somatic_indels_pass_final.txt
echo "#Variant must appear in at least one tumor samples after variant allele frequency and variant coverage filters" >> ./data/somatic_indels_pass_final.txt
echo "#Variant must have a coverage across entire cohort of tumor and germline bam files of greater than 150 base reads total (unique filter due to spots in genome of low coverage)" >> ./data/somatic_indels_pass_final.txt

# Header for filtered variants that did not pass 
echo -e "#Somatic indels that did not pass variant allele count filter\n" > ./data/somatic_indels_fail_allele_count.txt
echo -e "#Somatic indels that did not pass variant allele frequency filter\n" > ./data/somatic_indels_fail_vaf.txt
echo -e "#Somatic indels that did not pass tumor sample frequency threshold\n" > ./data/somatic_indels_fail_tumor_sample_freq.txt
echo -e "#Somatic indels that did not pass germline sample frequency threshold\n" > ./data/somatic_indels_fail_germline_sample_freq.txt

# Reading each unique somatic variant in while loop (input file at very end of loop)
while read line
do
chr=`echo $line | cut -d ' ' -f 1`
pos=`echo $line | cut -d ' ' -f 2`
ref=`echo $line | cut -d ' ' -f 3`
# Length of the reference allele
l_ref=`expr length $ref`
# check alt allele for comma, meaning 2 alt alleles. I customly checked these and one candidate of interest (FAN1) was a germline SNP
alt_int=`echo $line | cut -d ' ' -f 4`
if (echo $alt_int | grep -q ','); then
alt="Z"
else
alt=`echo $line | cut -d ' ' -f 4`
# Length of the alternative allele
l_alt=`expr length $alt`
fi
# Length of insertion according to samtools
# Arithmetic in POSIX shells is done with $ and double parentheses
# $(()) is preferable since it avoids a fork/execute for the expr command
l_st_ins="$(($l_alt - $l_ref))"
# Substring insertion from the right side, AKA the sequence of insertion in samtools mpileup
seq_st_ins=`echo -e "+$l_st_ins${alt:(-$l_st_ins)}"`
# Length of deletion according to samtools
l_st_del="$(($l_ref - $l_alt))"
l_st_del4seq=$((l_st_del+1))
seq_del=`seq -s "N" $l_st_del4seq | sed 's/[0-9]//g'`
seq_st_del=`echo -e "-$l_st_del$seq_del"`
# Variant sample
vs=`echo $line | cut -d ' ' -f 10`
# Variant with annotation
var=`echo $line`
# The gene symbol of gene associated with variant
symbol=`echo "$var" | cut -d " " -f 7`
# File name of variant
fv=`echo -e "$var" | sed 's/ /_/g'`
# Create file to count number of somatic samples containing variant, allows files to have a 0 count
touch "./data/filter/"$fv"_tumor.test"
# Create file to count number of germline samples containing variant, allows files to have a 0 count
touch "./data/filter/"$fv"_germline.test"
# Create file to count (2nd count) number of somatic samples containing variant, allows files to have a 0 count
touch ./data/filter/tumor_count.int
# Create file to count (2nd count) number of germline samples containing variant, allows files to have a 0 count
touch ./data/filter/germline_count.int
echo -e "\n$var"  >> ./data/filter/length.int
echo -e "REF: $ref"  >> ./data/filter/length.int
echo -e "ALT: $alt"  >> ./data/filter/length.int
echo -e "Insertion Length: $l_st_ins"  >> ./data/filter/length.int
echo -e "Insertion Sequence: $seq_st_ins"  >> ./data/filter/length.int
echo -e "Deletion Length: $l_st_del" >> ./data/filter/length.int
echo -e "Deletion Sequence: $seq_st_del"  >> ./data/filter/length.int
# Search each tumor bam for each somatic called insertions and deletions
# Create zero value for sum of coverage of tumor variant sites
tcovsum=0
# Create a zero value for the sum of variant allele frequency for each insertion
tvaf_ins_sum=0
# Create a zero value for the sum of variant allele frequency for each deletion
tvaf_del_sum=0
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/*-[1-9]*_S*_Bwa_RG_dedupped_realigned.bam`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
tmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# Tumor sample
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
tint=`echo $tmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
# Variant allele count
tvac_ins=`echo $tint | grep -o "$seq_st_ins" | wc -l`
tvac_del=`echo $tint | grep -o -- "$seq_st_del" | wc -l`
# Coverage
tcov=`echo -e "$tmpu" | cut -f4`
tcov="${tcov:-100000000000000000000000000000}"
# Add coverage to the sum of coverage
tcovsum=`expr $tcovsum + $tcov`
# Collect stats on variant counts of variant in 1 and all samples
tvaf_ins=`echo "scale=2;$tvac_ins/$tcov" | bc`
tvaf_ins="${tvaf_ins:-0}"
tvaf_ins_sum=`echo "scale=2;$tvaf_ins_sum + $tvaf_ins" | bc`
tvaf_del=`echo "scale=2;$tvac_del/$tcov" | bc`
tvaf_del="${tvaf_del:-0}"
tvaf_del_sum=`echo "scale=2;$tvaf_del_sum + $tvaf_del" | bc`

# Place somaticly verified variants in file so line number corresponds to number of tumor samples that contain varaint
if [ "$l_alt" -gt "1" ] && (echo $tint | grep -q "$seq_st_ins"); then
echo -e "$s\t$tmpu" >> "./data/filter/"$fv"_tumor.test"
elif [ "$l_ref" -gt "1" ] && (echo $tint | grep -q \\"$seq_st_del"); then
echo -e "$s\t$tmpu" >> "./data/filter/"$fv"_tumor.test"
fi
done
# Search each germline bam for each somatic called insertion and deletion
# Create zero value for sum of coverage of tumor variant sites
gcovsum=0
# Create a zero value for the sum of variant allele frequency for each insertion
gvaf_ins_sum=0
# Create a zero value for the sum of variant allele frequency for each deletion
gvaf_del_sum=0
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/*-[0]*_S*_Bwa_RG_dedupped_realigned.bam`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
gmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# Tumor sample
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
gint=`echo $gmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
# Variant allele count
gvac_ins=`echo $gint | grep -o "$seq_st_ins" | wc -l`
gvac_del=`echo $gint | grep -o -- "$seq_st_del" | wc -l`
# Coverage
gcov=`echo -e "$gmpu" | cut -f4`
gcov="${gcov:-100000000000000000000000000000}"
# Add coverage to the sum of coverage
gcovsum=`expr $gcovsum + $gcov`
# Collect stats on variant counts of variant in 1 and all samples
gvaf_ins=`echo "scale=2;$gvac_ins/$gcov" | bc`
gvaf_ins="${gvaf_ins:-0}"
gvaf_ins_sum=`echo "scale=2;$gvaf_ins_sum + $gvaf_ins" | bc`
gvaf_del=`echo "scale=2;$gvac_del/$gcov" | bc`
gvaf_del="${gvaf_del:-0}"
gvaf_del_sum=`echo "scale=2;$gvaf_del_sum + $gvaf_del" | bc`
# Place somatically called variants  that appear in germline samples in file so line number corresponds to number of germline samples that contain varaint
if [ "$l_alt" -gt "1" ] && (echo $gint | grep -q "$seq_st_ins"); then
echo -e "$s\t$gmpu" >> "./data/filter/"$fv"_germline.test"
elif [ "$l_ref" -gt "1" ] && (echo $gint | grep -q \\"$seq_st_del"); then
echo -e "$s\t$gmpu" >> "./data/filter/"$fv"_germline.test"
fi
done
# Create threshold value for variant allele count
vaf_ins_thres=`echo "scale=2;$gvaf_ins_sum" \* 2 | bc`
vaf_del_thres=`echo "scale=2;$gvaf_del_sum" \* 2 | bc`

# Add the sum of coverage on all tumor and germline samples. This will be used to filter out variants with inadequate
# coverage across cohort of samples
covsum=`expr $tcovsum + $gcovsum`
# Number of tumor samples with variant
tsnum=`wc -l "./data/filter/"$fv"_tumor.test" | cut -d ' ' -f1`
echo -e "$tsnum $fv" >> ./tumor_counts.txt
# Number of germline samples with variant
gsnum=`wc -l "./data/filter/"$fv"_germline.test" | cut -d ' ' -f1`
echo -e "$tsnum $fv" >> ./germline_counts.txt
# Print the header for variants in pass file
echo -e "\nVariant" | tee ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
echo -e "CHROM\tPOS\tREF\tALT\tMUT\tIMPACT\tSYMBOL\tENSEMBL\tENTREZ\tSAMPLE" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
echo "$var" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
echo "Tumor Samples Containing Variant" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
echo -e "SAMPLE\tCHROM\tPOS\tREF\tCOV\tBASES\tQUAL\tVAC\tVAF" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt

# Search each tumor bam for each somatic called variant
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/*-[1-9]*_S*_Bwa_RG_dedupped_realigned.bam`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
mpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# Tumor sample
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# CHROM, POS, READ_BASES
int=`echo $mpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
# FILTER VARIABLES
# Variant allele count
vac_ins=`echo $int | grep -o "$seq_st_ins" | wc -l`
vac_del=`echo $int | grep -o -- "$seq_st_del" | wc -l`
# Coverage assigned a value even if mpu variable is either empty or is undefined
cov=`echo -e "$mpu" | cut -f 4`
cov="${cov:-1000000}"
# Variant allele frequency even if mpu variable is either empty or is undefined
vaf_ins=`echo "scale=2;$vac_ins/$cov" | bc`
vaf_ins="${vaf_ins:-0}"
vaf_del=`echo "scale=2;$vac_del/$cov" | bc`
vaf_del="${vaf_del:-0}"
# Variant allele frequency threshold
vaf_thres=`echo "scale=2;1/20" | bc`

# Filters applying to insertions
if [ "$l_alt" -gt "1" ]; then
# Filter variants that contain less than 1 alternative allele
if [ "$vac_ins" -le "0" ]; then
echo -e "$s\t$mpu\t$vac_ins\t$vaf_ins" >> ./data/somatic_indels_fail_allele_count.txt
echo -e "Insertion Allele Count: $vac_ins" >> ./data/somatic_indels_fail_allele_count.txt
# Filter variants that contain a variant allele frequency less than variant allele frequency theshold (5%)
elif (( $(echo "$vaf_ins < $vaf_thres" | bc -l) )); then
echo -e "$s\t$mpu\t$vac_ins\t$vaf_ins" >> ./data/somatic_indels_fail_vaf.txt
echo -e "Variant Allele Frequency: $vaf_ins" >> ./data/somatic_indels_fail_vaf.txt
# Filter variant that do not show a variant allele in any tumor samples
# Redundent filter with current settings
elif [ "$tsnum" -le "0" ]; then
echo -e "$tsnum $fv" >> ./data/somatic_indels_fail_tumor_sample_freq.txt
# Filter variants that appear in greater than 8 germline samples
elif [ "$gsnum" -gt "8" ]; then
echo -e "$gsnum $fv" >> ./data/somatic_indels_fail_germline_sample_freq.txt
elif (( $(echo "$tvaf_ins_sum < $vaf_ins_thres" | bc -l) )); then
# Filter out variants without twice as many variants in somatic tissues compared to germline tissues across entire cohort
echo -e "$fv" >> ./data/filter/vaf.test
echo -e "\ntvaf_ins_sum: $tvaf_ins_sum" >> ./data/filter/vaf.test
echo "vaf_ins_thres: $vaf_ins_thres" >> ./data/filter/vaf.test
else
# Print the variables that passed all filters
echo -e "$s\t$mpu\t$vac_ins\t$vaf_ins" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
fi
fi
# Filters applying to deletions
if [ "$l_ref" -gt "1" ]; then
# Filter variants that contain less than 1 alternative allele
if [ $vac_del -le "0" ]; then
echo -e "$s\t$mpu\t$vac_del\t$vaf_del" >> ./data/somatic_indels_fail_allele_count.txt
echo -e "Deletion Allele Count: $vac_del" >> ./data/somatic_indels_fail_allele_count.txt
# Filter variants that contain a variant allele frequency less than variant allele frequency theshold (5%)
elif (( $(echo "$vaf_del < $vaf_thres" | bc -l) )); then
echo -e "$s\t$mpu\t$vac_del\t$vaf_del" >> ./data/somatic_indels_fail_vaf.txt
echo -e "Deletion Allele Frequency: $vaf_del" >> ./data/somatic_indels_fail_vaf.txt
# Filter variant that do not show a variant allele in any tumor samples
# Redundent filter with current settings
elif [ "$tsnum" -le "0" ]; then
echo -e "$tsnum $fv" >> ./data/somatic_indels_fail_tumor_sample_freq.txt
# Filter variants that appear in greater than 8 germline samples
elif [ "$gsnum" -gt "8" ]; then
echo -e "$gsnum $fv" >> ./data/somatic_indels_fail_germline_sample_freq.txt
elif (( $(echo "$tvaf_del_sum < $vaf_del_thres" | bc -l) )); then
# Filter out variants without twice as many variants in somatic tissues compared to germline tissues across entire cohort 
echo -e "$fv" >> ./data/filter/vaf.test
echo -e "\ntvaf_del_sum: $tvaf_del_sum" >> ./data/filter/vaf.test
echo "vaf_del_thres: $vaf_del_thres" >> ./data/filter/vaf.test
else
# Print the variables that passed all filters
echo -e "$s\t$mpu\t$vac_del\t$vaf_del" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
fi
fi
done

# Place holder to allow for counting of tumor samples with variants (2nd count)
if [ "$tsnum" -ge "1" ] && [ "$gsnum" -le "8" ]; then
echo -e "----------" | tee -a ./data/filter/tumor_count.int >> ./data/somatic_indels_pass.txt
fi

# Count the number of tumor samples that passed the test
# Print the variant as part of the header
echo "$var" >> indels_tumor_samples_pass.txt
# Print more header
echo "Tumor Samples That Passed Filters" >> indels_tumor_samples_pass.txt
# Count the number of lines between the strings "Tumor" and "---" in order to determine how many tumor samples for each varaint
tsnum2=`sed -n '/^Tumor/,/^---/p' ./data/filter/tumor_count.int | grep -v -e "^Tumor" -e "^---" -e "^SAMPLE" | wc -l`
# Capture the lines with the information about the variant calls
tscalls=`sed -n '/^Tumor/,/^---/p' ./data/filter/tumor_count.int | grep -v -e "^Tumor" -e "^---" -e "^SAMPLE"`
# Place a string at the beginning of the line for file manipulation later
tscalls_tsf=`echo -e "$tscalls" | sed "s/$/\t$tsnum2/g" | sed "s/^/LINE_T\t/"`
echo -e "$tsnum2\ntscalls_tsf" >> indels_tumor_samples_pass.txt

# Print the header and variants in FINAL pass file, (only if there is a tumor sample with called variant)
# For insertions
if [ "$l_alt" -gt "1" ] && [ "$tsnum" -ge "1" ] && [ "$gsnum" -le "8" ] && [ "$tsnum2" -ge "1" ] && [ "$covsum" -gt "150" ] && (( $(echo "$tvaf_ins_sum >= $vaf_ins_thres" | bc -l) )); then
echo -e "\nVariant" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_ID\tCHROM\tPOS\tREF\tALT\tMUT\tIMPACT\tSYMBOL\tENSEMBL\tENTREZ\tSAMPLE" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_V $var" >> ./data/somatic_indels_pass_final.txt
echo -e "Tumor Samples Containing Variant" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_ID\tSAMPLE\tCHROM\tPOS\tREF\tCOV\tBASES\tQUAL\tVAC\tVAF\tTSF" >> ./data/somatic_indels_pass_final.txt
echo -e "$tscalls_tsf" >> ./data/somatic_indels_pass_final.txt
fi
# For deletions
if [ "$l_ref" -gt "1" ] && [ "$tsnum" -ge "1" ] && [ "$gsnum" -le "8" ] && [ "$tsnum2" -ge "1" ] && [ "$covsum" -gt "150" ] && (( $(echo "$tvaf_del_sum >= $vaf_del_thres" | bc -l) )); then
echo -e "\nVariant" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_ID\tCHROM\tPOS\tREF\tALT\tMUT\tIMPACT\tSYMBOL\tENSEMBL\tENTREZ\tSAMPLE" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_V $var" >> ./data/somatic_indels_pass_final.txt
echo -e "Tumor Samples Containing Variant" >> ./data/somatic_indels_pass_final.txt
echo -e "LINE_ID\tSAMPLE\tCHROM\tPOS\tREF\tCOV\tBASES\tQUAL\tVAC\tVAF\tTSF" >> ./data/somatic_indels_pass_final.txt
echo -e "$tscalls_tsf" >> ./data/somatic_indels_pass_final.txt
fi
done <./data/all_nonsyn_indels_to_read.txt

# Make a copy of the pass file for further file manipulation
cp ./data/somatic_indels_pass_final.txt ./data/somatic_indels_pass_final2.txt
# Capture the lines that start with these 2 strings and change the delimiter from space to tab
grep -e "^LINE_V" -e "^LINE_T" ./data/somatic_indels_pass_final2.txt | \
sed 's/ /\t/g' > ./data/somatic_indels_pass_final3.txt

# Create an intermediate file called layout to mimic final file
echo -e "#CHROM\tPOS\tREF\tALT\tCOV\tBASES\tVAC\tVAF\tMUT\tIMPACT\tSYMBOL\tENSEMBL\tENTREZ\tSAMPLES" > ./data/layout.txt
# Read each line of the file with variants
while read line
do
# Apply the proper values from each line (String ID's at the beginning of each line to allow easy sorting)
start=`echo -e "$line" | cut -f1`
if [ "$start" = "LINE_V" ]; then
v_line=`echo -e "$line"`
chr=`echo "$v_line" | cut -f2`
pos=`echo "$v_line" | cut -f3`
ref=`echo "$v_line" | cut -f4`
alt=`echo "$v_line" | cut -f5`
mut=`echo "$v_line" | cut -f6`
impact=`echo "$v_line" | cut -f7`
symbol=`echo "$v_line" | cut -f8`
ensembl=`echo "$v_line" | cut -f9`
entrez=`echo "$v_line" | cut -f10`
elif [ "$start" = "LINE_T" ]; then
t_line=`echo -e "$line"`
samples=`echo -e "$t_line" | cut -f2`
cov=`echo -e "$t_line" | cut -f6`
bases=`echo -e "$t_line" | cut -f7`
vac=`echo -e "$t_line" | cut -f9`
vaf=`echo -e "$t_line" | cut -f10`
else
:
fi
if [ "$start" = "LINE_T" ]; then
echo -e "$chr\t$pos\t$ref\t$alt\t$cov\t$bases\t$vac\t$vaf\t$mut\t$impact\t$symbol\t$ensembl\t$entrez\t$samples" >> ./data/layout.txt
fi
done <./data/somatic_indels_pass_final3.txt

#Sort for unique values
(grep "^#" ./data/layout.txt; grep -v "^#" ./data/layout.txt | sort | uniq) > ./data/layout_uniq.txt

# Perform additional filtering step with variant allele frequency
while read line
do
if (echo $line | grep -q "^#CHROM"); then
echo "$line" > ./data/layout_uniq2.txt
else
alt=`echo "$line" | cut -f4`
l_alt=`expr length $alt`
ref=`echo "$line" | cut -f3`
l_ref=`expr length $ref`
tbird_sample=`echo "$line" | cut -f14`
tbird=`echo "$line" | cut -f14 | cut -c 1-3`
chr=`echo "$line" | cut -f1`
pos=`echo "$line" | cut -f2`
cov=`echo "$line" | cut -f5`
# For insertions
if [ "$l_alt" -gt "1" ]; then
# Length of insertion according to samtools
# Arithmetic in POSIX shells is done with $ and double parentheses
# $(()) is preferable since it avoids a fork/execute for the expr command
l_st_ins="$(($l_alt - $l_ref))"
# Substring insertion from the right side, AKA the sequence of insertion in samtools mpileup
seq_st_ins=`echo -e "+$l_st_ins${alt:(-$l_st_ins)}"`
# Length of deletion according to samtools
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/$tbird-[0]*_S*_Bwa_RG_dedupped_realigned.bam`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
gmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# Germline sample
gsample=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
gint=`echo $gmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
gcall=`echo -e "$gsample\t$gmpu"`
# Variant allele count
gvac_ins=`echo $gint | grep -o "$seq_st_ins" | wc -l`
gvaf_ins=`echo "scale=2;$gvac_ins/$cov" | bc`
gvaf_ins="${gvaf_ins:-0}"
gvaf_ins_thres=`echo "scale=2;$gvaf_ins" \* 2 | bc`
done
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/$tbird_sample"_Bwa_RG_dedupped_realigned.bam"`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
tmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# tumor sample
tsample=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
int=`echo $tmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
tcall=`echo -e "$tsample\t$tmpu"`
# Variant allele count
tvac_ins=`echo $int | grep -o "$seq_st_ins" | wc -l`
tvaf_ins=`echo "scale=2;$tvac_ins/$cov" | bc`
tvaf_ins="${tvaf_ins:-0}"
done
if (( $(echo "$tvaf_ins >= $gvaf_ins_thres" | bc -l) )) && [ "$l_alt" -gt "1" ]; then
echo -e "$line" >> ./data/layout_uniq2.txt
#echo "tvaf_ins: $tvaf_ins" >> ./data/layout_uniq2.txt
#echo "gvaf_ins: $gvaf_ins" >> ./data/layout_uniq2.txt
#echo "gvaf_ins_thres: $gvaf_ins_thres" >> ./data/layout_uniq2.txt
fi
fi
# For deletions
if [ "$l_ref" -gt "1" ]; then
# Length of deletion according to samtools
# Arithmetic in POSIX shells is done with $ and double parentheses
# $(()) is preferable since it avoids a fork/execute for the expr command
# Length of deletion according to samtools
l_st_del="$(($l_ref - $l_alt))"
l_st_del4seq=$((l_st_del+1))
seq_del=`seq -s "N" $l_st_del4seq | sed 's/[0-9]//g'`
seq_st_del=`echo -e "-$l_st_del$seq_del"`
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/$tbird-[0]*_S*_Bwa_RG_dedupped_realigned.bam`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
gmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# Germline sample
gsample=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
gint=`echo $gmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
gcall=`echo -e "$gsample\t$gmpu"`
# Variant allele count
gvac_del=`echo $gint | grep -o -- "$seq_st_del" | wc -l`
gvaf_del=`echo "scale=2;$gvac_del/$cov" | bc`
gvaf_del="${gvaf_del:-0}"
gvaf_del_thres=`echo "scale=2;$gvaf_del" \* 2 | bc`
done
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/$tbird_sample"_Bwa_RG_dedupped_realigned.bam"`
do
# CORE VARIABLES
# CHROM, POS, REF, COV, BASES, QUAL
tmpu=`samtools mpileup -r $chr:$pos-$pos $bam`
# tumor sample
tsample=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# Samtools mpileup output: CHROM, POS, READ_BASES
tint=`echo $tmpu | cut -d ' ' -f 1,2,5 | tr '[:lower:]' '[:upper:]'`
tcall=`echo -e "$tsample\t$tmpu"`
# Variant allele count
tvac_del=`echo $tint | grep -o -- "$seq_st_del" | wc -l`
tvaf_del=`echo "scale=2;$tvac_del/$cov" | bc`
tvaf_del="${tvaf_del:-0}"
done
if (( $(echo "$tvaf_del >= $gvaf_del_thres" | bc -l) )) && [ "$l_ref" -gt "1" ]; then
echo "$line" >> ./data/layout_uniq2.txt
#echo "tvaf_del: $tvaf_del" >> ./data/layout_uniq2.txt
#echo "gvaf_del: $gvaf_del" >> ./data/layout_uniq2.txt
#echo "gvaf_del_thres: $gvaf_del_thres" >> ./data/layout_uniq2.txt
fi
fi
fi
done <./data/layout_uniq.txt

# Sort for the unique variants
grep -v "^#" ./data/layout_uniq2.txt | cut -f1,2,3,4 | sort | uniq > ./data/uniq_vars.int

# Grep the layout file for uniq variants and collect the appropriate stats
rm ./results/somatic_indels_final.int
echo -e "#CHROM\tPOS\tREF\tALT\tMUT\tIMPACT\tSYMBOL\tENSEMBL\tENTREZ\tTSN\tSAMPLE\tVAC\tVAF" > ./results/somatic_indels_final.txt
while read line
do
# Set value for variatn allele count
vac=`grep -e "$line" ./data/layout_uniq2.txt | cut -f7`
vacs=`echo -e "$vac" | tr '\n' '|' | sed 's/|$//'`
vaf=`grep -e "$line" ./data/layout_uniq2.txt | cut -f8`
vafs=`echo -e "$vaf" | tr '\n' '|' | sed 's/|$//'`
sample=`grep -e "$line" ./data/layout_uniq2.txt | cut -f14 | awk '!seen[$0]++'`
samples=`echo -e "$sample" | tr '\n' '|' | sed 's/|$//'`
tsn=`echo -e "$samples" | tr '|' '\n' | wc -l`
part1=`grep -e "$line" ./data/layout_uniq2.txt | cut -f1-4 | awk '!seen[$0]++'`
p1n=`echo $part1 | wc -l`
part2=`grep -e "$line" ./data/layout_uniq2.txt | cut -f9-13 | awk '!seen[$0]++' | head -n $p1n`
echo -e "$part1\t$part2\t$tsn\t$samples\t$vacs\t$vafs" >> ./results/somatic_indels_final.int
done <./data/uniq_vars.int
sort -k7,7 ./results/somatic_indels_final.int >> ./results/somatic_indels_final.txt


# Perform a quick check to see which genes are mutated by both somatic SNVs and Indels
for gene in `cut -f7 ./results/somatic_snvs_final.txt | grep -v -e "SYMBOL" -e "NA" -e "LOC" | sort | uniq`
do
grep "$gene" ./results/somatic_indels_final.txt
done

# Perform samtools manual analyses for each potential variant
# This will be used in conjunction with manual analysis of each gene to determine further filtering techniques
rm ./results/double_check.txt
while read line
do
echo -e "\n$line" >> ./results/double_check.txt
chr=`echo -e "$line" | cut -f1`
pos=`echo -e "$line" | cut -f2`

echo -e "\nTUMOR SAMPLES" >> ./results/double_check.txt
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/*-[1-9]*_S*_Bwa_RG_dedupped_realigned.bam`
do
mpu=`samtools mpileup -r $chr:$pos-$pos $bam`
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
echo -e "$s\t$mpu" >> ./results/double_check.txt
done

echo -e "\nGERMLINE SAMPLES" >> ./results/double_check.txt
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/*-[0]*_S*_Bwa_RG_dedupped_realigned.bam`
do
mpu=`samtools mpileup -r $chr:$pos-$pos $bam`
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
echo -e "$s\t$mpu" >> ./results/double_check.txt
done

echo -e "\nPARENTAL SAMPLES" >> ./results/double_check.txt
for bam in `ls -1 /home/proj/MDW_genomics/xu/final_bam/00268[3-4]_Line-[6-7]_Bwa_RG_dedupped_realigned.bam`
do
mpu=`samtools mpileup -r $chr:$pos-$pos $bam`
s=`echo "$bam" | xargs -i basename {} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
echo -e "$s\t$mpu" >> ./results/double_check.txt
done
done <./results/somatic_indels_final.txt

