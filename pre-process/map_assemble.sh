#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 6 
#$ -M
########################

# Parameters
declare -a all_sam=('D18-6449' 'D18-6450' 'D18-6451' 'D18-6452' 'D18-6453' 'D18-6454'
		    'D18-6455' 'D18-6456' 'D18-6458' 'D18-6459' 'D18-6460' 'D18-6461'
		    'D18-6463' 'D18-6464' 'D18-6465' 'D18-6466' 'D18-6468' 'D18-6469'
		    'D18-6470' 'D18-6472');

sample_group=all_sam

# Set up environment
activate=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin/activate;
hisat2=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/envs/hisat2/;
export PATH=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin:$PATH
export PATH=/bin:$PATH
source $activate $hisat2;

echo 'Virtual environment active...'
cd /net/bmc-pub7/data0/essigmannlab/jobs/hisat2

# Run alignments using HISAT2
echo 'Running HISAT2...'
date

for s in ${all_sam[@]}; do
  echo 'Aligning sample' $s'...'
  hisat2 --dta -x ref/mm10_mod -U trimmed/$s/180711Ess\_$s.fastq -S $s.hisat2.sam
  samtools view -bS $s.hisat2.sam > $s.hisat2.unsorted.bam
  samtools sort -o $s.hisat2.bam $s.hisat2.unsorted.bam
done

echo 'Alignments complete.'
date

# Assemble transcripts using HISAT2
echo 'Assembling transcripts for each sample...'
date

for s in ${all_sam[@]}; do
  echo 'Assembling sample' $s'...'
  stringtie -G ref/mm10.gtf -A $s\_gene\_abund.tab -o $s.gtf -l $s $s.hisat2.bam
done
echo 'Stringtie (re-)assembly complete.'
date

echo 'Merging transcripts from all samples...'
stringtie --merge -G ref/mm10.gtf -o $sample_group\_merged.gtf $sample_group\_mergelist.txt

# Compare transcripts with reference
echo 'Comparing transcripts with reference...'
gffcompare -r ref/mm10.gtf -G -o merged $sample_group\_merged.gtf

# Estimate transcript abundances and table counts
echo 'Estimating abundances for differential expression analysis...'
mkdir ballgown
for s in ${all_sam[@]}; do
  echo 'Running abundance estimates for sample' $s'...'
  mkdir ballgown/$s
  stringtie -e -B -G $sample_group\_merged.gtf -o ballgown/$s/$s.gtf $s.hisat2.bam
done

echo 'Abundance estimates complete.'
date
