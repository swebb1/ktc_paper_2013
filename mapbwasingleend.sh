
echo "bwa aln reads..."

bwa aln -t 15 /homes/genomes/s.cerevisiae/sacCer3/bwa_indexes/sacCer3 $1 > $2.sai

echo "bwa samse..."

bwa samse /homes/genomes/s.cerevisiae/sacCer3/bwa_indexes/sacCer3 $2.sai $1 > $2_se.sam

echo "Mapping complete, filtering for uniquely aligned reads..."

echo "Sam to Bam..."

samtools view -hSb $2_se.sam -o $2_se.bam

echo "Bam sort..."

samtools sort $2_se.bam $2_se_sort

echo "Index bam..."

samtools index $2_se_sort.bam

echo "Flag stats..."

samtools flagstat $2_se_sort.bam > $2_se_sort.flagstat

echo "Filter bam..."

samtools view -hb -q 20 -F 4 $2_se_sort.bam -o $2_se_sort_filter.bam

echo "Index Filtered bam..."

samtools index $2_se_sort_filter.bam

echo "Remove duplicates..."

java -jar /storage/software/picard-tools-1.54/MarkDuplicates.jar I=$2_se_sort_filter.bam O=$2_se_sort_filter_NR.bam M=$2.dupmetrics REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE

echo "Index NR bam..."

samtools index $2_se_sort_filter_NR.bam

echo "Flag stats..."

samtools flagstat $2_se_sort_filter.bam > $2_se_sort_filter.flagstat

echo "Flag stats..."

samtools flagstat $2_se_sort_filter_NR.bam > $2_se_sort_filter_NR.flagstat

rm $2_se.sam $2_se.bam
