# rna-seq
Pipeline for RNA-seq scripts used by the Essigmann Lab.

## Part I: Pre-Process, Align, and Assemble

### Setup
1. Create environment from `yml` file: `conda env create -f rnaseq_env.yml`
2. Activate environment: `source activate rnaseq`

### Prepare FASTA reference
1. Download genome from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/): `chromFa.tar.gz`
2. Unzip FASTA: `tar -xvzf chromFa.tar.gz`
3. Remove mitochondrial chromosome and other noncanonical chromosomes (`chr#_#########_random`) from directory
4. Compile chromosomal FASTA files to single file: `cat *.fa > mm10.fa`
5. If necessary, modify FASTA file to match naming convention for GTF file: `sed -i 's/chr//g' mm10.fa`
6. Index reference: `hisat2-build -f mm10.fa mm10`

### Prepare GTF reference transcriptome
1. Download GTF from [Ensembl](https://bit.ly/2xPCJYz): `Mus_musculus.GRCm38.93.gtf.gz`
2. Unzip GTF: `tar -xzvf Mus_musculus.GRCm38.93.gtf.gz`
3. (Optional) Rename file: `mv Mus_musculus.GRCm38.93.gtf mm10.gtf`
4. Format known splice junctions to format used by HISAT2: `hisat2_extract_splice_sites.py mm10.gtf > mm10.gtf.ss`

### Trim raw RNA-seq reads
1. Trim adapter sequences and ends: `trimmomatic-0.38.jar SE $seq.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:30 LEADING:30 TRAILING:30 MINLEN:25`

### Align RNA-seq reads to reference genome using HISAT2
1. Map to whole genome, accounting for known splice sites: `hisat2 --dta -x ref/mm10 --known-splicesite-infile ref/mm10.gtf.ss -U $trimmed.fastq -S $sample.hisat2.sam`
2. Convert to BAM: `samtools view -bS $sample.hisat2.sam > $sample.hisat2.unsorted.bam`
3. Sort BAM file: `samtools sort -o $sample.hisat2.bam $sample.hisat2.unsorted.bam`

### Assemble and quantify expressed genes and transcripts with StringTie
1. Estimate abundances for differential expression analysis: `stringtie -e -B -G ref/mm10.gtf -A $sample\_abund.tab -o $directory/$sample/$sample.gtf $sample.hisat2.bam`
   * _Note:_ This is considered StringTie's "alternate" workflow, relying on a well-annotated reference; it will not search for novel isoforms. Suggested by the StringTie creator [here](https://github.com/gpertea/stringtie/issues/170).
   * _Historical note:_ Originally the pipeline used StringTie's recommended workflow, but identifying gene names caused trouble as many `gene_id` values were given `MSTRG` assignments. StringTie author made note of it [here](https://github.com/gpertea/stringtie/issues/179).

### Prepare StringTie outputs for differential expression analysis
1. Download Python script ([prepDE.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)) provided by StringTie developers
2. Run script to extract read count information from StringTie outputs: `python2.7 prepDE.py -l 40 [-i $directory]`
   * _Note:_ Without the `-i` parameter, this assumes default directory structure created by StringTie, with a `ballgown` folder in the working directory. In our case, use the `-i` parameter to denote the directory where outputs are contained.
   * _Note:_ The script requires a Python version between 2.7 and 3.
   * _Note:_ The `-l` parameter takes in average read length. While this doesn't affect relative transcript levels, it will impact your absolute values! The default parameter for `-l` is `75`.
   
## Part II: Normalization and Differential Gene Analysis

### Setup
1. R environment must have the following Bioconductor packages installed: DESeq2, org.Mm.eg.db, and biomaRt.
   * If not installed, install Bioconductor: `source("https://bioconductor.org/biocLite.R")`
   * Install packages: `biocLite("DESeq2"); biocLite("org.Mm.eg.db"); biocLite("biomaRt")`
   * Call each library in the R environment with the `library()` function.
2. Read in expression gdata calculated by StringTie: `pheno_data <- read.table("180711Ess_phenodata.csv", header=TRUE, sep=","); count_data <- read.table("gene_count_matrix.csv", header=TRUE, sep=",", row.names=1)`
3. Retrieve gene IDs and gene names by creating a `biomaRt` instance: `ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")`
4. Create key-value pairs for Ensembl IDs to gene names: `gene_names <- getBM(c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id", values=rownames(count_data), mart=ensembl); colnames(gene_names) <- c("ensembl_id", "gene_name"); gene_set <- setNames(as.list(gene_names$gene_name), as.list(gene_names$ensembl_id))`
5. Read in all data to a `DESeqDataSet`: `dds <- DESeqDataSetFromMatrix(countData=count_data, colData=pheno_data, design~1)`
6. Filter to remove low-count genes: `keep.rows <- rowSums(counts(dds)) >= 10; dds <- dds[keep.rows,]`

### Conduct and output differential gene analysis
1. Use in-house R function `dea.group` to query between either sexes or treatment groups. For sex differences, the function requires the DESeq data structure and experimental treatment group to query; TCPOBOP addition changes require the DESeq data structure,the base group (i.e. non-TCPOBOP treatment), and the sex to query.
   * Conduct differential gene analysis on DMSO-treated mice, by sex: `dds.sex_dmsona <- dea.group(dds, "sex", "DMSO")`
   * Conduct differential gene analysis on DMSO-treated males, by TCPOBOP addition: `dds.tc_dmsom <- dea.group(dds, "TCPOBOP", "DMSO", "M")`

## References
* [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic)
* [_Nature Protocols_ paper on HISAT2 and StringTie](https://www.nature.com/articles/nprot.2016.095)
* [HISAT2 & StringTie manual](https://ccb.jhu.edu/software/hisat2/manual.shtml)
