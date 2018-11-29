#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 6
#$ -M klkim@mit.edu
##############################

# Parameters
declare -a all_sam=('D18-6449' 'D18-6450' 'D18-6451' 'D18-6452' 'D18-6453' 'D18-6454'
		    'D18-6455' 'D18-6456' 'D18-6458' 'D18-6459' 'D18-6460' 'D18-6461'
		    'D18-6463' 'D18-6464' 'D18-6465' 'D18-6466' 'D18-6468' 'D18-6469'
		    'D18-6470' 'D18-6472');

sample_group=all_sam
runID=180711Ess

# Set up environment
activate=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin/activate;
rnaseq=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/envs/rnaseq/;
export PATH=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin:$PATH
export PATH=/bin:$PATH
source $activate $rnaseq;

echo 'Virtual environment active...'
cd /net/bmc-pub7/data0/essigmannlab/jobs/rnaseq

# Trim adapter sequences and low-quality ends
data_path=/net/bmc-pub7/data0/essigmannlab/data/$runID
trim_path=/net/bmc-pub7/data0/essigmannlab/dependencies/Trimmomatic-0.38

mkdir 01
for s in ${all_sam[@]}; do
  echo 'Trimming sample' $s'...'
  java -jar $trim_path/trimmomatic-0.38.jar SE $data_path/$s/$runID\_$s\_NA_sequence.fastq 01/$runID\_$s.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:30 LEADING:30 TRAILING:30 MINLEN:25
done

echo 'Trimming complete.'
date

# Run alignments using HISAT2
echo 'Running HISAT2...'
date

mkdir 02
for s in ${all_sam[@]}; do
  echo 'Aligning sample' $s'...'
  hisat2 --dta -x ref/mm10 -U 01/$s/$runID\_$s.fastq -S 02/$s.hisat2.sam
  samtools view -bS 02/$s.hisat2.sam > 02/$s.hisat2.unsorted.bam
  samtools sort -o 02/$s.hisat2.bam 02/$s.hisat2.unsorted.bam
done

echo 'Alignments complete.'
date

# Estimate transcript abundances and table counts
echo 'Estimating abundances for differential expression analysis...'
date

mkdir 03
for s in ${all_sam[@]}; do
  echo 'Running abundance estimates for sample' $s'...'
  mkdir 03/$s
  stringtie -e -B -G ref/mm10.gtf -A 03/$s\_abund.tab -o 03/$s/$s.gtf 02/$s.hisat2.bam
done

echo 'Abundance estimates complete.'
date

# Create gene count matrix
echo 'Creating gene count matrix...'
date

python2.7 scripts/prepDE.py -l 40 -i 03

echo 'Gene count matrix created.'
date

