# rna-seq
Pipeline for RNA-seq scripts used by the Essigmann Lab.

## Part I: Pre-Process, Align, and Assemble

### Setup
1. Create environment from `yml` file: `conda env create -f rnaseq_env.yml`
2. Activate environment: `source activate rnaseq`

### Prepare FASTA reference and GTF
1. Download genome from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/): `chromFA.tar.gz`
2. Download GTF from [Ensembl](https://bit.ly/2xPCJYz): `Mus_musculus.GRCm38.93.gtf.gz`
3. Unzip FASTA: `tar -xvzf chromFa.tar.gz`
4. Remove mitochondrial chromosome and other noncanonical chromosomes (`chr#_#########_random`) from directory
5. Compile chromosomal FASTA files to single file: `cat *.fa > mm10.fa`
6. If necessary, modify FASTA file to match naming convention for GTF file: `sed -i 's/chr//g' mm10.fa`
7. Index reference: `hisat2-build -f mm10.fa mm10`

### Trim raw RNA-seq reads
1. Trim adapter sequences and ends: `trimmomatic-0.38.jar SE $seq.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:30 LEADING:30 TRAILING:30 MINLEN:25`

## References
* [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic)
* [Nature Protocols paper on HISAT2 and StringTie](https://ccb.jhu.edu/software/hisat2/manual.shtml)
* [HISAT2 manual](https://www.nature.com/articles/nprot.2016.095)
