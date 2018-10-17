# rna-seq
Pipeline for RNA-seq scripts used by the Essigmann Lab.

## Part I: Pre-Process, Align, and Assemble

### Setup
1. Create environment from `yml` file: `conda env create -f rnaseq_env.yml`
2. Activate environment: `source activate rnaseq`

### Prepare FASTA reference and GTF
1. Download genome from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/): `chromFA.tar.gz`
2. Unzip FASTA: `tar -xvzf chromFa.tar.gz`
3. Remove mitochondrial chromosome and other noncanonical chromosomes (`chr#_#########_random`) from directory
4. Compile chromosomal FASTA files to single file: `cat *.fa > mm10.fa`
5. If necessary, modify FASTA file to match naming convention for GTF file: `sed -i 's/chr//g' mm10.fa`
6. Index reference: `hisat2-build -f mm10.fa mm10`

### Prepare GTF transcriptome reference
1. Download GTF from [Ensembl](https://bit.ly/2xPCJYz): `Mus_musculus.GRCm38.93.gtf.gz`
2. Unzip GTF: `tar -xzvf Mus_musculus.GRCm38.93.gtf.gz`
3. (Optional) Rename file: `mv Mus_musculus.GRCm38.93.gtf mm10.gtf`
4. Format known splice junctions to format used by HISAT2: `hisat2_extract_splice_sites.py mm10.gtf > mm10.gtf.ss`

### Trim raw RNA-seq reads
1. Trim adapter sequences and ends: `trimmomatic-0.38.jar SE $seq.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:30 LEADING:30 TRAILING:30 MINLEN:25`

### Align RNA-seq reads to reference genome using HISAT2
1. Map to whole genome, accounting for known splice sites: `hisat2 --dta -x ref/m10 --known-splicesite-infile ref/mm10.gtf.ss -U $trimmed.fastq -S $sample.hisat2.sam`
2. Convert to BAM: `samtools view -bS $sample.hisat2.sam > $sample.hisat2.unsorted.bam`
3. Sort BAM file: `samtools sort -o $sample.hisat2.bam $sample.hisat2.unsorted.bam`

### Assemble and quantify expressed genes and transcripts with StringTie
1. Assemble transcripts: `stringtie -G ref/mm10.gtf [-A $sample\_gene\_abund.tab] -o $sample.gtf -l $sample $sample.hisat2.bam`
2. Merge transcripts from all samples: `stringtie --merge -G ref/mm10.gtf -G -o merged merged\_mergelist.txt`
3. Compare sequenced transcripts with the reference: `gffcompare -r ref/mm10.gtf -G -o merged merged.gtf`
4. Estimate abundances for differential expression analysis: `stringtie -e -B -G merged.gtf -o ballgown/$sample/$sample.gtf $sample.hisat2.bam`

### Prepare StringTie outputs for differential expression analysis
1. Download Python script ([prepDE.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)) provided by StringTie developers
2. Run script to extract read count information from StringTie outputs: `python prepDE.py`
   * Note: This assumes default directory structure created by StringTie, with a `ballgown` folder in the working directory.
   * Note: The script requires a Python version between 2.7 and 3.

## References
* [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic)
* [Nature Protocols paper on HISAT2 and StringTie](https://ccb.jhu.edu/software/hisat2/manual.shtml)
* [HISAT2 manual](https://www.nature.com/articles/nprot.2016.095)
