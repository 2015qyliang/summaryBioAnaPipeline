published paper: [**2017**-Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life](https://www.nature.com/articles/s41564-017-0012-7)

# Methods

## 01-Recovery of cultivation-independent genomes

Metadata for metagenomes in the [Sequence Read Archive (SRA)](https://academic.oup.com/nar/article/39/suppl_1/D19/2505848) at the National Center for Biotechnology Information (NCBI) were obtained from the [SRAdb (within R)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-19).

Only metagenomes submitted to the SRA **before 31 December 2015** were **considered** with a predominant **focus on environmental** and non-human gastrointestinal samples (for example, rumen, guinea pigs and baboon faeces; Supplementary Table 1).

Metagenomes from studies where MAGs had previously been recovered were excluded if the UBA MAGs did not provide appreciable improvements in genome quality or phylogenetic diversity.

Each of the 1,550 metagenomes were processed independently, with all SRA Runs within an SRA Experiment (that is, sequences from a single biological sample) being co-assembled using the **CLC de novo assembler v.4.4.1 (CLCBio**) [ref-1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4701411/) [ref-2](https://web.uri.edu/gsc/files/Genomics-Workbench-Guide.pdf) [ref-3](http://resources.qiagenbioinformatics.com/tutorials/De_novo_assembly_paired_data.pdf).

Assembly was restricted to contigs ≥500 bp and the word size, bubble size and paired-end insertion size determined by the assembly software.

Assembly statistics are reported for contigs ≥2,000 bp (Supplementary Table 1).

Reads were mapped to contigs with [**BWA v.0.7.12-r1039**](https://academic.oup.com/bioinformatics/article/25/14/1754/225615) (*this big*) using the **BWA-MEM algorithm** with default parameters and the mean coverage of contigs obtained using the ‘coverage’ command of [CheckM v.1.0.6](https://genome.cshlp.org/content/25/7/1043).

Genomes were independently recovered from each SRA Experiment using [MetaBAT v.0.26.3](https://peerj.com/articles/1165/) *under all five preset parameter settings (that is, verysensitive, sensitive, specific, veryspecific, superspecific)*.

The completeness and contamination of the genomes recovered under each MetaBAT preset were estimated using CheckM using lineage-specific markers genes and default parameters.

For each SRA Experiment, only genomes recovered with the MetaBAT preset resulting in the largest number of bins with an estimated completeness >70% and contamination <5% were considered for further refinement and validation.

## 02-Merging of compatible bins

Automated binning methods can produce multiple bins from the same microbial population.

The ‘merge’ method of **CheckM v.1.0.6** was used to identify pairs of bins where the completeness increased by ≥10% and the contamination increased by ≤1% when **merged into a single bin**.

Bins meeting these criteria were grouped into a single bin if the **mean GC of the bins were within 3%**, the **mean coverage of the bins had an absolute percentage difference ≤25%**, and the bins had identical taxonomic classifications as determined by their placement in the reference genome tree used by CheckM.

This set of criteria was used to avoid producing chimaeric bins.

## 03-Filtering scaffolds with **divergent** genomic properties

Scaffolds with genomic features deviating substantially from the mean GC, tetranucleotide signature, or coverage of a bin were identified with the ‘outliers’ method of [**RefineM v.0.0.14**](https://github.com/dparks1134/RefineM) using default parameters.

This removes all scaffolds with a GC or tetranucleotide distance outside the **98th percentile** of the expected distributions of these genomic features, as determined empirically over a set of 5,656 trusted reference genomes [CheckM](https://genome.cshlp.org/content/25/7/1043),[Methane metabolism in the archaeal phylum Bathyarchaeota revealed by genome-centric metagenomics](http://science.sciencemag.org/content/350/6259/434).

 Scaffolds were also removed if their mean coverage had an absolute percentage **difference ≥50%** when compared to the mean coverage of the bin.

## 04-Filtering scaffolds with **incongruent taxonomic classification**

Each gene within a bin was assigned a taxonomic classification through homology search using [BLASTP v.2.2.30+](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421) against a custom database of 12,321 genomes from [**RefSeq/GenBan release 75**](https://academic.oup.com/nar/article/42/D1/D553/1066302).

This database was constructed from RefSeq and GenBank genomes consisting of **≤300 contigs**, having an **N50 ≥20 kb** and **containing ≤10 kb of ambiguous base pairs**.

A genome was only included in the database if it was estimated to be ≥90% complete, ≤10% contaminated and had an overall quality ≥50 (defined as completeness − 5 × contamination).

**Quality estimates** were determined with [**CheckM**](http://ecogenomics.github.io/CheckM/manual/checkm_manual.pdf) using the [**lineage-specific workflow**](https://github.com/Ecogenomics/CheckM/wiki/Workflows) [*biostars*](https://www.biostars.org/p/195935/) and default parameters.

Genomes meeting this set of requirements were dereplicated to remove genomes from the same named species with an **amino-acid identity (AAI) ≥99.5%**. [example AAI](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132499/)

**AAI** values were calculated with [**CompareM v.0.0.13**](https://github.com/dparks1134/CompareM) and dereplication performed in a greedy fashion with a preference towards type strains and genomes annotated as complete at NCBI.

Genes were** assigned the taxonomic classification** of their ‘top’ hit or designated as unclassified if the gene had no identified homologue with an **E-value ≤1e−2**, a percent **sequence identity ≥30%** and a percent **alignment length ≥50%**.

Scaffolds with **incongruent taxonomic classifications** were removed from each bin.

The consensus classification of a bin at each taxonomic **rank** was determined by identifying the taxon that occurred at the highest frequency across all classified genes or designated as unclassified if no taxon was represented by ≥50% of the classified genes.

Scaffolds where **≥50%** of the classified genes at each rank agreed with the consensus classification of the bin were designated as ‘trusted’, and a taxon was considered to be ‘common’ if it comprised ≥5% of the classified genes across the set of trusted scaffolds.

A scaffold was considered to be taxonomically incongruent and removed from a bin if the following three conditions were met:
- (1) it contained ≥5 classified genes and ≥25% of all genes on the scaffold were classified;
- (2) ≤10% of the classified genes were contained in the set of common taxa at each classified rank;
- (3) >50% of classified genes were assigned to the same taxon at each classified rank.

Taxonomic classification of genes and identification of scaffolds with divergent taxonomic classifications were performed with the **‘taxon_profile’** and **‘taxon_filter’** methods of [**RefineM v.0.0.14**](https://github.com/dparks1134/RefineM), respectively.

## 05-Filtering scaffolds with incongruent **16S rRNA genes**

Scaffolds were removed from a bin if they contained a complete or partial 16S rRNA gene **≥600 bp** with a taxonomic classification incongruent with the taxonomic identity of the bin.

[**BLASTN**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421) was used to assign 16S rRNA genes the taxonomy of its closest homologue within a database comprising the 10,769 16S genes identified within the 12,321 reference genomes discussed in the previous section.

The sequence identity to the closest **homologue** was used to determine the set of ranks that should be examined for congruency.

Specifically, previously reported median percent identities values were used to establish [**conservative thresholds for the taxonomic ranks to consider: genus ≥98.7%, family ≥96.4%, order ≥92.25%, class ≥89.2%, phylum ≥86.35% and domain ≥83.68%**](https://www.nature.com/articles/nrmicro3330) (**2014** Uniting the classification of cultured and uncultured bacteria and archaea using 16S rRNA gene sequences).

The taxon at each rank was then compared to the taxonomic classification of the genes across all scaffolds in the bin and designated as incongruent if the taxon was assigned to ≤10% of classified genes.

This methodology is implemented in [**RefineM v.0.0.14**](https://github.com/dparks1134/RefineM).

## 06-Selection of refined genomes

Of the 64,295 bins produced by **MetaBAT**, only the 7,903 genomes with an estimated quality ≥50 (defined as completeness − 5 × contamination), scaffolds resulting in an **N50 of ≥10 kb**, **containing <100 kb ambiguous bases** and consisting of **<1,000 contigs** and **<500 scaffolds** were considered to be of sufficient quality for further exploration and deposition in public repositories.

We adopted the quality criteria of completeness − 5 × contamination as it provides a good signal (completeness) to noise (contamination) ratio, where higher levels of contamination are only permissible when the genome is largely complete.

These genomes have been deposited as assemblies in NCBI’s TPA:Assembly database along with alignment files indicating the mapping of SRA reads to UBA genomes.

## 07-Comparison of UBA genomes to complete conspecific strains

The 3,438 near-complete (≥90% complete; ≤5% contamination) UBA genomes were compared to complete isolate genomes in RefSeq release 76.

Of these, 207 of the UBA genomes were determined to be conspecific strains of complete isolate genomes based on an **average nucleotide identity (ANI)** and **alignment fraction (AF)** above [96.5% and 60%, respectively](https://academic.oup.com/nar/article/43/14/6761/2903001) -- *Microbial species delineation using whole genome sequences*. ANI and AF values were determined using [**ANI Calculator v.1**]((https://academic.oup.com/nar/article/43/14/6761/2903001)).

The genome size of the UBA genomes was adjusted to account for its estimated completeness and contamination: **adjusted genome size** = **(genome size)/(completeness + contamination)**.

**Homologues** between UBA genomes and their conspecific counterparts were determined by inferring genes with [Prodigal v.2.6.3](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) and establishing sequence **similarity** with **BLASTP** v.2.2.30+.

A UBA protein was **considered homologous** to an isolate protein if it was the top hit among all isolate proteins, had an **E-value of ≤1e−10**, a percent **identity of ≥70%** and an **alignment length spanning ≥70%** of the isolate protein.

## 08-Proteins used to infer genome trees

Bacterial and archaeal genome trees were inferred from the concatenation of 120 (Supplementary Table 6) and 122 (Supplementary Table 7) phylogenetically informative proteins, respectively.

These proteins were identified as being present in ≥90% of bacterial or archaeal genomes and, when present, single-copy in ≥95% of genomes.

**Protein-coding regions** were identified using **Prodigal v.2.6.3** (with default parameters, but with *Ns treated as masked sequences*), translation tables determined using a [**coding density heuristic**](https://genome.cshlp.org/content/25/7/1043) (CheckM), and the **ubiquity of genes** determined across genomes from NCBI’s RefSeq release 73 annotated with the [Pfam v.27](https://academic.oup.com/nar/article/42/D1/D222/1062431) and [TIGRFAMs v.15.0 databases](https://academic.oup.com/nar/article/31/1/371/2401559).

Only genomes composed of** ≤200 contigs**, with an **N50 of ≥20 kb** and with CheckM **completeness and contamination estimates of ≥95% and ≤5%**, respectively, were considered.

Phylogenetically informative proteins were determined by **filtering ubiquitous proteins** whose gene trees had poor congruence with a set of subsampled concatenated genome trees.

Specifically, the initial set of 188 bacterial (187 archaeal) proteins were randomly subsampled to 132 genes (~70%) and concatenated to infer a subsampled genome tree.

**Gene subsampling** was independently performed 100 times to establish well-supported splits, which we define as any split occurring in >80% of the subsampled trees and with ≥1% of taxa contained in both bipartitions induced by the split.

*other links:[The effects of subsampling gene trees on coalescent methods applied to ancient divergences](https://www.sciencedirect.com/science/article/pii/S1055790315003929)-[bioconductor-subSeq](https://rdrr.io/bioc/subSeq/man/subsample.html)-[subSeq](https://academic.oup.com/bioinformatics/article/30/23/3424/207085)-[Subsample high coverage single cells](https://jdblischak.github.io/singleCellSeq/analysis/subsample-high-coverage-lcl.html)*

The congruence between a gene tree and the subsampled genome tree was measured as the fraction of well-supported split lengths compatible with the gene tree, a measure we call the **‘normalized compatible split length’**.

Genes with a normalized compatible split length of ≤50% were removed, as this poor congruence may indicate the presence of lateral gene transfer events.
http://darwin.uvigo.es/download/prottest_manual.pdf
**Proteins** were aligned to Pfam and TIGRfam HMMs using [**HMMER v.3.1b1**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195) with default parameters and trees were inferred with [**FastTree v.2.1.7**](https://academic.oup.com/mbe/article/26/7/1641/1128976) under the **WAG+GAMMA** models.

*other links:[ProtTest: Selection of best-fit models of protein evolution](http://darwin.uvigo.es/download/prottest_manual.pdf)--[RAxML hands-on session](https://sco.h-its.org/exelixis/web/software/raxml/hands_on.html)--[Trends in Amino Acid Substitution Models](https://www.frontiersin.org/articles/10.3389/fgene.2015.00319/full)*

Trees were also inferred from two ribosomal protein sets:
- (1) **16 ribosomal proteins** (Supplementary Table 4) that form a **syntenic block** [ref--A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648) [10--Unusual biology across a group comprising more than 15% of domain Bacteria](https://www.nature.com/articles/nature14486)
- (2) **23 ribosomal proteins** (Supplementary Table 5) previously used for tree inference and tested for [**lateral gene transfer** -- Insights into the phylogeny and coding potential of microbial dark matter](https://www.nature.com/articles/nature12352).

## 09-Inference of genome trees

Genome trees were inferred across a dereplicated set of UBA and RefSeq/GenBank release 76 (May 2016; includes 727 single cell and 1,811 MAGs) genomes.

All 5,192 RefSeq genomes annotated at NCBI as ‘reference’ or ‘representative’ were retained, except for a low-quality subset of 294 genomes that did not meet **our ‘trusted’ genome criteria**: composed of ≤300 contigs, N50 ≥20 kb, CheckM completeness and contamination estimates of ≥90% and ≤10%, respectively.

This set of 4,898 genomes was augmented with an additional 3,324 RefSeq genomes to retain at least two genomes per species where possible.

Preference was given to genomes annotated at NCBI as being a type strain and/or ‘complete’ and restricted to genomes meeting the ‘trusted’ genome filtering criteria.

An additional 551 RefSeq genomes currently without a species designation at NCBI, but passing the genome quality filtering, were also added to this initial set of seed genomes.

UBA, GenBank and remaining RefSeq genomes meeting the ‘trusted’ genome criteria were compared to these 8,773 seed genomes.

**Genomes with an AAI of ≥99.5% to a seed genome**, as calculated over the 120 bacterial or 122 archaeal marker genes used for phylogenetic inference, were clustered with the seed genome and do not appear as separate genomes in the genome trees.

This **cutoff correlated with the proposed 96.5% ANI threshold** for [defining bacterial and archaeal species](https://academic.oup.com/nar/article/43/14/6761/2903001) (Supplementary Fig. 6).

Trusted genomes with an **AAI <99.5% **were added to the seed set.

All remaining genomes, regardless of quality, were compared to this final seed set using the **same AAI clustering criteria of 99.5%**.

**Seed and unclustered genomes** with an estimated genomes quality ≥50 (defined as completeness − 5 × contamination) were used to create an initial multiple sequence alignment, with the exception of the 797 CPR genomes10, which were retained regardless of their estimated quality.

Proteins were **identified and aligned** using HMMER v.3.1b1 and the resulting **alignment trimmed** to remove columns represented by <50% of taxa or without a common amino acid in ≥25% of taxa.

Genomes with **amino acids in <40%** of aligned columns (20% for the lenient archaeal trees) were **removed** from consideration.

The 120 concatenated bacterial protein set consisted of 34,796 aligned columns after trimming and was inferred over 19,198 genomes.

The 122 concatenated archaeal protein set contained 28,025 trimmed columns and spanned 1,012 genomes when using standard filtering criteria and 27,942 columns spanning 1,070 genomes when using lenient filtering.

Trees were inferred with **FastTree v.2.1.7 under the WAG+GAMMA models** and support values determined using **100 non-parametric bootstrap** replicates.

## 10-Inference of 16S rRNA trees

Bacterial and archaeal trees were inferred from 16S rRNA genes >600 bp and >1,200 bp within UBA and RefSeq/GenBank release 76 genomes, respectively.

The 16S rRNA genes were identified using **HMMER** and domain-specific SSU/LSU HMM models as implemented in the ‘ssu-finder’ method of CheckM. These genes were aligned with [**ssu-align v.0.1**](https://academic.oup.com/bioinformatics/article/25/10/1335/270663) and trailing or leading columns represented by ≤70% of taxa trimmed, which resulted in bacterial and archaeal alignments of 1,421 and 1,378 bp, respectively.

Trees were inferred with **FastTree** v.2.1.7 under the **GTR+GAMMA models** and support values determined using 100 non-parametric bootstrap replicates.

## 11-Similarity of 16S rRNA genes

The percent identity between 16S rRNA genes was calculated from the multiple sequence alignments used to infer the domain-specific 16S rRNA gene trees.

The **‘dist.seqs’** command of [**mothur v.1.30.2**](https://aem.asm.org/content/75/23/7537) was used to calculate percent identity.

Default parameters were used, **except** that gaps at the end of sequences were ignored (countends = F) in order to accommodate partial 16S rRNA sequences.

*Inter-phylum (inter-class) 16S rRNA percent identity values were determined by identifying the most similar sequence to each sequence within a phylum across all sequences from different phyla (classes)*.

## 12-Assessment of phylogenetic and taxonomic diversity

**Phylogenetic diversity** (*total branch length spanned by a set of taxa*) and **gain** (additional branch length contributed by a set of taxa) were calculated using [**GenomeTreeTk v.0.0.23**](https://github.com/dparks1134/GenomeTreeTk) and verified with [**ARB v.6.0.2**](https://academic.oup.com/nar/article/32/4/1363/1038459).

**Taxonomic diversity** and the percentage of lineages of equal evolutionary distance unique to the UBA genomes were determined using the mean branch length to extant [taxa criterion -- A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648).

**Lineages of equal evolutionary distance** were related to the distribution of NCBI taxa as defined on 19 May 2016 and used to construct the phylum-level lineage view (Supplementary Fig. 7) by evaluating the number of groups formed at **mean branch length values of 0.5 to 1.1 with a step size of 0.025**.

A value of 0.85 was selected as it most closely matched the number of bacterial phyla when excluding the CPR.

In agreement with previous analyses, we used this criterion to explore the taxonomic structure of phylogenetic trees and not to explicitly establish taxonomic status.

## 13-Genomic similarity

The **AAI and shared gene content** between genomes were determined with **CompareM v.0.0.21**(https://github.com/dparks1134/CompareM) using default parameters (**homologues** defined by an E-value ≤0.001, a percent identity ≥30% and an alignment length ≥70%).

[**CompareM**](https://github.com/dparks1134/CompareM) reports **shared gene content** relative to the genome with the fewest identified genes in order to accommodate incomplete genomes.

Inter-phylum and inter-class AAI and shared gene content values were determined by sampling up to 50 near-complete (completeness ≥90%, contamination ≤5%, N50 >20 kb, total contigs ≤200) RefSeq release 76 genomes from each named lineage, taking care to sample evenly between named species.

The **AAI score**, defined as the **sum of the AAI and shared gene content**, was used to determine the most similar genome to each query genome.

## 14-Genomic and assembly properties

Genomic and assembly properties (for example, GC, N50, coverage) were determined using **CheckM**. Transfer RNAs were identified with [tRNAscan-SE v.1.3.1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC146525/) using either the bacterial or archaeal tRNA model and default parameters.
