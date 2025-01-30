name='pATM_seb_all'
size="2000"

#bgzip -d wgEncodeAwgTfbsUtaGm12878CtcfUniPk.narrowPeak.gz
Rscript --vanilla get_fasta.R  pATMDIvAseb_cleaned_q0.001_peaks.narrowPeak ${name} ${size}
cut -f 1-3 hg19_blacklisted/${name}_hg19_blacklisted.bed > hg19_blacklisted/$name.cut.bed

sed -i 's/"//g' hg19_blacklisted/$name.cut.bed
sed -i "s/$/\tpATM/" hg19_blacklisted/$name.cut.bed
echo "sed ok"
sort -k1V -k2n -k3n hg19_blacklisted/$name.cut.bed > hg19_blacklisted/sorted_$name.bed
echo "sort ok"
bgzip -c hg19_blacklisted/sorted_$name.bed > hg19_blacklisted/sorted_$name.bed.gz
echo "bgzip ok"
tabix -p bed hg19_blacklisted/sorted_$name.bed.gz
echo "tabix ok"
cut -f 4 hg19_blacklisted/sorted_$name.bed | sort -u > hg19_blacklisted/distinct_features.txt

cut -f 1-3 hg19_blacklisted_400bestAsiSI_masked/${name}hg19_blacklisted_400bestAsiSI_masked_${size}.bed > hg19_blacklisted_400bestAsiSI_masked/$name.cut.bed

sed -i 's/"//g'  hg19_blacklisted_400bestAsiSI_masked/$name.cut.bed

sed -i "s/$/\tpATM/" hg19_blacklisted_400bestAsiSI_masked/$name.cut.bed
echo "sed ok"
sort -k1V -k2n -k3n hg19_blacklisted_400bestAsiSI_masked/$name.cut.bed > hg19_blacklisted_400bestAsiSI_masked/sorted_$name.bed
echo "sort ok"
bgzip -c hg19_blacklisted_400bestAsiSI_masked/sorted_$name.bed > hg19_blacklisted_400bestAsiSI_masked/sorted_$name.bed.gz
echo "bgzip ok"
tabix -p bed hg19_blacklisted_400bestAsiSI_masked/sorted_$name.bed.gz
echo "tabix ok"
cut -f 4 hg19_blacklisted_400bestAsiSI_masked/sorted_$name.bed | sort -u > hg19_blacklisted_400bestAsiSI_masked/distinct_features.txt
