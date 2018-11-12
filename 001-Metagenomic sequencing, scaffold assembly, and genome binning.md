# Metagenomic sequencing, scaffold assembly, and genome binning

Published paper: [Genomic expansion of magnetotactic bacteria reveals an early common origin of magnetotaxis with lineage-specific evolution](https://www.nature.com/articles/s41396-018-0098-9)

## Sequencing

Metagenomic DNA was extracted and amplified from magnetically enriched MTB cells as previously described [Origin of microbial biomineralization and magnetotaxis during the Archean](http://www.pnas.org/content/114/9/2171). Shotgun sequencing of metagenomic DNA from each location was performed with an Illumina HiSeq 2000 using the pair-end 2×125 reads with a 600-bp insert size or using the Illumina HiSeq 4000 with the pair-end strategy of 150-bp reads with an average 270-bp insert size (Beijing Genomics Institute, Beijing, China)

##  Scaffold assembly

Illumina reads were trimmed to remove the adapter sequences and low-quality bases, and were assembled using metaSPAdes [metaSPAdes: a new versatile metagenomic assembler](https://genome.cshlp.org/content/early/2017/03/15/gr.213959.116) with the following parameters (--only-assembler -k 31,41,51,61,71,81,91,101,111,121)

## Genome binning

Assembled scaffolds ≥2500 bp were binned separately using MetaBAT v0.26.1 [MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities](https://peerj.com/articles/1165/) and MyCC [Accurate binning of metagenomic contigs via automated clustering sequences using information of genomic signatures and marker genes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4828714/). Results of two binning methods for each sample were combined and a non-redundant set of bins was chosen. 

The acquired genomes were curated manually with two approaches: 
- 1) using the CheckM [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/) “outliers” command to identify scaffolds from bins that appear to be outliers in either GC, tetranucleotide, or coding density space relative to the expected distribution of these genomic statistics;
- 2) using BLASTn or BLASTx to identify potential contaminant contigs based on their top BLAST hits.

## Assess

The quality and accuracy of the acquired genomes were assessed using CheckM [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/) based on the taxonomic-specific workflow (domain Bacteria) and QUAST [QUAST: quality assessment tool for genome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/)

## Annotation

Genomes were annotated using Prokka v1.11 [Prokka: rapid prokaryotic genome annotation](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) with manual curation.  The average amino-acid identity (AAI) was calculated using enveomics [The enveomics collection: a toolbox for specialized analyses of microbial genomes and metagenomes](https://peerj.com/preprints/1900/).



