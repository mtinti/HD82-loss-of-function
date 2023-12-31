#!/bin/bash
#$ -l local_free=100G,m_mem_free=20G,h_rt=43200
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -N s_1_Ald4
#$ -o /cluster/majf_lab/mtinti/Anna_DH/F21FTSEUHT0831_TRYxdpqR/Ald4/s_1_Ald4_$JOB_ID

set -e
set -u 
set -o pipefail

experiment='Ald4'
genome='tb927_10'
base_fastq='Ald4_'

echo 'copy indata files'
#copy experiment folders in TMPDIR
mkdir -p $TMPDIR'/genomes/'$genome'/' && cp -avr  'genomes/'$genome'/' $TMPDIR'/genomes/' 
mkdir -p $TMPDIR'/'$experiment'/data/' && cp -avr  $experiment'/data/' $TMPDIR'/'$experiment'/'

#direct link to genome files
path_genome_index=$TMPDIR'/genomes/'$genome'/'$genome
path_trascriptome_index=$TMPDIR'/trascriptome/'$genome'/'$genome

echo 'folder $TMPDIR'
ls -l $TMPDIR

echo 'make out dirs'
#folders to store output files
path_out=$TMPDIR'/'$experiment'/'
mkdir $path_out'fastqc/'
mkdir $path_out'fastp'
mkdir $path_out'qualimap_bam'
mkdir $path_out'qualimap_rna'
mkdir $path_out'featureCounts'
mkdir $path_out'bowtie'
mkdir $path_out'markduplicates' 
mkdir $path_out'flagstat' 
mkdir $path_out'stats' 
mkdir $path_out'alfred'
mkdir $path_out'unmap'
mkdir $path_out'fastqs'

ls $path_out -l
##qc of fastq 

#direct link to fastq files
fastq_1=$path_out'data/'$base_fastq'1.fq.gz'
fastq_2=$path_out'data/'$base_fastq'2.fq.gz'


#test to check fastq file consistency
for f in $(ls $TMPDIR'/'$experiment'/data/'*gz); do echo -ne "$f:\t" && c=$(zcat $f|wc -l) && echo -ne "$c\t" && echo $c%4|bc; done


echo 'find size of fastq file'
inFileLength=$(echo $(zcat $fastq_1 | wc -l)/4|bc)
echo 'inFileLength all: '$inFileLength

echo 'run fastqc' 
fastqc  $fastq_1 -o $path_out'fastqc/'
fastqc  $fastq_2 -o $path_out'fastqc/'

#remove problematic reads
echo 'run fastp' 
out_1=$path_out'data/'$base_fastq'f1.fastq.gz'
out_2=$path_out'data/'$base_fastq'f2.fastq.gz'
fastp -i $fastq_1 -I $fastq_2 -o $out_1 -O $out_2 -h $path_out'fastp/fastp_'$base_fastq'.html' -j $path_out'fastp/'$base_fastq'_fastp.json'


base_fastq=$experiment'_'

echo 'run bowtie2'
(bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $out_1 -2 $out_2) \
2>$path_out'bowtie/'$experiment'.log' | samtools view -bSu -q 0 | \
samtools sort -@ 8 -o $path_out$base_fastq'sorted.bam'

samtools index $path_out$base_fastq'sorted.bam'


echo 'run picard MarkDuplicates' 
picard MarkDuplicates \
I=$path_out$base_fastq'sorted.bam' \
O=$path_out$base_fastq'sortedd.bam' \
M=$path_out'markduplicates/marked_dup_metrics.txt' \
REMOVE_DUPLICATES=false \
ASSUME_SORT_ORDER=coordinate

mv $path_out$base_fastq'sortedd.bam' $path_out$base_fastq'sorted.bam'

samtools sort -@ 8 -o $path_out$base_fastq'sorted.bam' $path_out$base_fastq'sorted.bam'
samtools index $path_out$base_fastq'sorted.bam'

echo 'run alfred'
alfred qc -r $path_genome_index.fa -j $path_out'alfred/'$base_fastq'qc.json.gz' -o $path_out'alfred/'$base_fastq'qc.tsv.gz' $path_out$base_fastq'sorted.bam'


echo 'run samtools flagstat'
samtools flagstat -@ 8 $path_out$base_fastq'sorted.bam' > $path_out'flagstat/'$experiment'.txt'
echo 'run samtools stats'
samtools stats -@ 8 $path_out$base_fastq'sorted.bam' > $path_out'stats/'$experiment'.txt'


unset DISPLAY
echo 'run qualimap bamqc' 
qualimap bamqc --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam' \
-outdir $path_out'qualimap_bam/'$experiment'/' \
-outfile $experiment'.bamqc.html' -outformat 'HTML'

echo 'run qualimap rnaseq' 
qualimap rnaseq --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam'  \
-gtf $path_genome_index'.gtf' \
-outdir $path_out'qualimap_rna/'$experiment'/' \
-pe -outfile $experiment'.bamqc.html' -outformat 'HTML'


echo 'run genome coverage'
bedtools genomecov -ibam $path_out$base_fastq'sorted.bam' \
-bg -pc > $path_out$base_fastq'sorted_pc_bg.bed'
gzip $path_out$base_fastq'sorted_pc_bg.bed'


bamCoverage --numberOfProcessors 8 --binSize 100 \
-b $path_out$base_fastq'sorted.bam' -o $path_out$base_fastq'sorted_pc_bg.bw'


echo 'run featureCounts' 
#echo 'run 11'
#B BothEndsMapped
#C exclude chimeric (reads map on different chr)
#M count multi mapping
#O count overlapping
#T cores
#-t use the feature
#-g summarize by
#-a is a gtf file
featureCounts -p -B -C -M -O -T 8 -t transcript -g gene_id -a $path_genome_index'.gtf' \
-o $path_out'counts.txt' $path_out$base_fastq'sorted.bam'

cp $path_out'counts.txt.summary' $path_out'featureCounts/'


echo 'extract unmapped'
samtools view -bu -f 12 -F 256  $path_out$base_fastq'sorted.bam' \
| samtools sort -n -o $path_out$base_fastq'sortedname_unmap.bam'

echo 'extract fastq files'
bamToFastq -i $path_out$base_fastq'sortedname_unmap.bam' \
-fq $path_out'unmap/'$base_fastq'.sorted_unmap_all_1.fastq' -fq2 $path_out'unmap/'$base_fastq'.sorted_unmap_all_2.fastq'


echo 'setup fastq screen: copy genome'
mkdir -p $TMPDIR'/FastQ_Screen_Genomes/' && cp -avr  'FastQ_Screen_Genomes/' $TMPDIR

echo 'setup fastq screen: alter config file'
sed -i "s@tmp_folder@$TMPDIR@g" $TMPDIR'/FastQ_Screen_Genomes/fastq_screen2.conf'

echo 'run fastq screen'
fastq_screen --nohits $path_out'unmap/'$base_fastq'.sorted_unmap_all_1.fastq' $path_out'unmap/'$base_fastq'.sorted_unmap_all_2.fastq' \
--conf $TMPDIR'/FastQ_Screen_Genomes/fastq_screen2.conf' --outdir $path_out'fastqs' --aligner bowtie2

echo 'gzip unmapped / remove unmapped'
rm -r $path_out'unmap/'

#gzip $path_out'unmap/'$base_fastq'.sorted_unmap_all_1.fastq'
#gzip $path_out'unmap/'$base_fastq'.sorted_unmap_all_2.fastq'
#B_pol_3_sortedname_unmap.bam

rm $path_out'fastqs/'$base_fastq'.sorted_unmap_all_1.tagged.fastq'
rm $path_out'fastqs/'$base_fastq'.sorted_unmap_all_2.tagged.fastq'

rm $path_out'fastqs/'*.fastq


#rm $path_out$base_fastq'sorted.bam'
#rm $path_out$base_fastq'sorted.bam.bai'
rm $path_out$base_fastq'sortedname_unmap.bam'
#rm $path_out$base_fastq'sorted_pc_bg.bed.gz'

echo 'copy results' 
rm -r $TMPDIR'/'$experiment'/data'
rm -r $TMPDIR'/FastQ_Screen_Genomes'

#rm $path_out$base_fastq*.bam
#rm $path_out$base_fastq*.bai
mkdir -p $experiment'/res2' && cp -avr $TMPDIR'/'$experiment $experiment'/res2'










