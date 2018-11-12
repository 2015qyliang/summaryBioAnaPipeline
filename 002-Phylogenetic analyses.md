# Phylogenetic analyses

Published paper: [Genomic expansion of magnetotactic bacteria reveals an early common origin of magnetotaxis with lineage-specific evolution](https://www.nature.com/articles/s41396-018-0098-9)

## Construct genomic tree

The maximum likelihood phylogeny of genomes was constructed using [RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies](https://academic.oup.com/bioinformatics/article/30/9/1312/238053) (-m PROTGAMMAVT -f a -x 12345 -k -p 12345 -N 100) with a concatenated alignment of the conserved **ubiquitous proteins** identified with [PhyloPhlAn is a new method for improved phylogenetic and taxonomic placement of microbes](https://www.nature.com/articles/ncomms3304)  

The VT+G model was used as determined by [ProtTest 3: fast selection of best-fit models of protein evolution](https://academic.oup.com/bioinformatics/article/27/8/1164/227935) 

The genomic tree was rooted with genomes from the domain Archaea (Methanobrevibacter ruminantium and Methanobrevibacter smithii). Confidence in phylogenetic results was assessed using the rapid bootstrap algorithm of RAxML with 100 replicates [A rapid bootstrap algorithm for the RAxML Web servers](https://academic.oup.com/sysbio/article/57/5/758/1618491). Bootstrap convergence test was conducted using RAxML (-I autoMRE). 

In order to further identify whether the Magnetococcales order represent a novel class in the Proteobacteria phylum, we additionally constructed a maximum likelihood phylogenomic tree with a concatenated **amino-acid sequence** alignment (6988 amino-acid positions) of 43 **lineage-specific marker genes** from 15 Magnetococcales genomes and up to 248 Proteobacteria genomes generated using the “tree” command in  [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/) 

The genomic tree was constructed using RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)  under the LG+I+G model of evolution.

## Homologous sequence analysis

**Homologous sequences** of magnetosome proteins MamA, -B, -E, -K, -M, and -Q were identified within the refseq_protein database using PSI-BLAST searches (BLOSUM62 scoring matrix, E-value < 1e-05, with exclusion of published MTB genomes) with each magnetosome protein from Magnetospirillum gryphiswaldense MSR-1, Desulfovibrio magneticus RS-1, and “Candidatus Magnetoglobus multicellularis” as query sequences. 

The hits were combined and clustered using CD-HIT [Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences](https://academic.oup.com/bioinformatics/article/22/13/1658/194225) with a sequence similarity cutoff of 0.8. The complete **amino-acid sequences** of magnetosome proteins MamA, -B, -E, -K, -M, and -Q from all available MTB genomes and their non-MTB homologs were aligned by MUSCLE: multiple sequence alignment with high accuracy and high throughput](https://academic.oup.com/nar/article/32/5/1792/2380623) algorithms using [MEGA](https://www.megasoftware.net/). 

Phylogenetic trees were then generated using the maximum likelihood method of [ FastTree 2—approximately maximum-likelihood trees for large alignments](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) with default settings. 

Multiple alignments of MamE, -M, and -Q were concatenated and a phylogenetic tree was constructed using RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phyloDgenies](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)  with the LG+I+G model as determined by [ProtTest 3: fast selection of best-fit models of protein evolution](https://academic.oup.com/bioinformatics/article/27/8/1164/227935). Confidence in a phylogenetic tree was assessed using 100 bootstrap replicates. 

## Trees visualiziton

Trees were visualized using [ FastTree 2—approximately maximum-likelihood trees for large alignments](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) and iTOL [ Interactive tree of life (iTOL)v3: an online tool for the display and annotation of phylogenetic and other trees](https://academic.oup.com/nar/article/44/W1/W242/2499315).




