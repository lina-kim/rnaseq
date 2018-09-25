# rna-seq
Pipeline for RNA-seq scripts used by the Essigmann Lab.

## Part I: Pre-Process, Align, and Assemble

### Setup
1. Create environment from `yml` file: `conda env create -f rnaseq_env.yml`
2. Activate environment: `source activate rnaseq`

### Prepare FASTA reference and GTF
1. Download genome from UCSC Genome Browser: `chromFA.tar.gz`
2. Download GTF from Ensembl: `Mus_musculus.GRCm38.93.gtf.gz`
3. Unzip FASTA: `tar -xvzf chromFa.tar.gz`
4. Remove mitochondrial chromosome and other noncanonical chromosomes (`chr#_#########_random`) from directory
