Published paper: [Nitrogen-fixing populations of Planctomycetes and Proteobacteria are abundant in surface ocean metagenomes](https://www.nature.com/articles/s41564-018-0176-9)

* heterotrophic bacterial diazotrophs (HBDs)*
* nitrogen fixation and nitrogenase-like sequences

# Metagenomic co-assemblies, gene calling and binning

We organized the data set into 12 ‘metagenomic sets’ based on the geographic coordinates of metagenomes [(Supplementary Table 1)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0176-9/MediaObjects/41564_2018_176_MOESM3_ESM.xlsx). We co-assembled reads from each metagenomic set using [MEGAHIT53 v1.0.3](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) [Github](https://github.com/voutcn/megahit), with a **minimum scaffold length of 1 kbp**, and simplified the scaffold header names in the resulting assembly outputs using [anvi’o v2.3.0](https://merenlab.org/software/anvio).
For each metagenomic set, we then binned scaffolds >2.5 kbp (>5 kbp for the Southern Ocean) following the workflow outlined in [anvi’o](https://peerj.com/articles/1319/).
- (1) anvi’o was used to profile the scaffolds using [Prodigal v2.6.3](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) with default parameters to identify genes [(Supplementary Table 2)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0176-9/MediaObjects/41564_2018_176_MOESM4_ESM.xlsx), and [HMMER55](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195) to identify genes matching to archaeal [56](https://www.nature.com/articles/nature12352) and bacterial [57](https://www.nature.com/articles/nmeth.3103),[58](http://www.pnas.org/content/110/14/5540),[59](https://www.nature.com/articles/ismej2011189),[60](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022099) single-copy core gene collections;
- (2) [Centrifuge](https://genome.cshlp.org/content/26/12/1721) was used with NCBI's NT database to infer the taxonomy of genes (as described in https://merenlab.org/2016/06/18/importing-taxonomy);
- (3) short reads were mapped from the metagenomic set to the scaffolds using [Bowtie2 v2.0.5](https://www.nature.com/articles/nmeth.1923) and the recruited reads stored as BAM files using [samtools](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) ;
- (4) anvi’o was used to profile each BAM file to estimate the coverage and detection statistics of each scaffold, and to combine mapping profiles into a merged profile database for each metagenomic set.

We then clustered scaffolds with the automatic binning algorithm [CONCOCT](https://www.nature.com/articles/nmeth.3103) by constraining the number of clusters per metagenomic set to 100 to minimize the ‘fragmentation error’ (when multiple clusters describe one population), with the exception of the Southern Ocean (25 clusters) and the Pacific Ocean southeast (150 clusters) metagenomic sets. Finally, we manually binned each CONCOCT cluster (n = 1,175) using the anvi’o interactive interface. [Supplementary Table 10](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0176-9/MediaObjects/41564_2018_176_MOESM12_ESM.xlsx) reports the genomic features (including completion and redundancy values) of the characterized bins.

# Identification and curation of MAGs

We defined all bins with >70% completeness or >2 Mbp in length as MAGs [(Supplementary Table 2)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0176-9/MediaObjects/41564_2018_176_MOESM4_ESM.xlsx). We then individually refined each MAG as outlined in [ref](https://peerj.com/articles/1839/), and renamed scaffolds they contained accordingly to their MAG ID to ensure that the names of all scaffolds in MAGs we characterized from the 12 metagenomic sets were unique

# Taxonomic and functional inference of MAGs

We used [CheckM](https://genome.cshlp.org/content/25/7/1043) to infer the taxonomy of MAGs based on the proximity of 43 single-copy gene markers within a reference genomic tree.
We also used Centrifuge, RAST and manual BLAST searches of single-copy core genes against the NCBI's non-redundant database to manually refine the CheckM taxonomic inferences, especially regarding the archaeal and eukaryotic MAGs.
We also used the occurrence of bacterial single-copy core genes to identify MAGs affiliated to the Candidate Phyla Radiation (as described in https://merenlab.org/2016/04/17/predicting-CPR-Genomes/).
We used KEGG to identify functions and pathways in MAGs.
We also used RAST to identify functions in 15 MAGs that contained the complete set of nitrogen fixation genes (originally identified from the KEGG pathways).
We used [Gephi](https://gephi.org/publications/gephi-bastian-feb09.pdf) v0.8.2 to generate a functional network using the Force Atlas 2 algorithm to connect MAGs and RAST functions. Node sizes were correlated to the number of edges they contained, which resulted in larger nodes for MAGs compared to functions.

# Characterization of a non-redundant database of MAGs

We concatenated all scaffolds from the genomic database of MAGs into a single FASTA file and used **Bowtie2 and samtools** to recruit and store reads from the 93 metagenomes.
We used **anvi’o** to determine the coverage values, detection and relative distribution of MAGs and individual genes across metagenomes.
The **Pearson correlation coefficient** of each pair of MAGs was calculated based on their relative distribution across the 93 metagenomes using the function ‘cor’ in R.
Finally, [**NUCmer**](https://academic.oup.com/nar/article/30/11/2478/1024948) was used to determine the **average nucleotide identity (ANI)** of each pair of MAGs affiliated to the same phylum for improved performance (the Proteobacteria MAGs were further split at the class level).
MAGs were considered redundant when their ANI reached 99% (minimum alignment of >75% of the smaller genome in each comparison) and the **Pearson correlation coefficient was above 0.9**.
We then selected a single MAG to represent a group of redundant MAGs based on the largest **‘completion minus redundancy’** value from single-copy core genes for Archaea and Bacteria, or longer genomic length for Eukarya.
This analysis provided a non-redundant genomic database of MAGs. We performed a final mapping of all metagenomes to calculate the mean coverage and detection of these MAGs

# **Statistical analyses**
[STAMP](https://academic.oup.com/bioinformatics/article/26/6/715/245265) and Welch’s test were used to identify non-redundant MAGs that were significantly enriched in the Pacific Ocean compared to all the other regions combined. Supplementary Table 3 reports the P values for each MAG.

# World maps

We used the ggplot2 package for R to visualize the metagenomic sets and relative distribution of MAGs in the world map.

# Phylogenomic analysis of MAGs

We used [PhyloSift v1.0.1](https://peerj.com/articles/243/) with default parameters to infer associations between MAGs in a phylogenomic context. Briefly
- (1) identifies a set of **37 marker gene** families in each genome,
- (2) concatenates the alignment of each marker gene family across genomes,
- (3) computes a phylogenomic tree from the concatenated alignment using [FastTree v2.1](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490).
We rooted the phylogenomic tree to the phylum Planctomycetes with [FigTree v1.4.3](http://tree.bio.ed.ac.uk/software/figtree/), and used [anvi’o](https://merenlab.org/software/anvio) to visualize it with additional data layers.
