published paper: [A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648)

# Methods

A data set comprehensively covering the **three domains of life** was **generated** using publicly available genomes **from the Joint Genome Institute's IMG-M database** (img.jgi.doe.gov), a previously developed data set of [eukaryotic genome information](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3768317/), previously published genomes derived from metagenomic data sets [7--Unusual biology across a group comprising more than 15% of domain Bacteria](https://www.ncbi.nlm.nih.gov/pubmed/26083755),[8--Genomic Expansion of Domain Archaea Highlights Roles for Organisms from New Phyla in Anaerobic Carbon Cycling](https://doi.org/10.1016/j.cub.2015.01.014),[31--Critical biogeochemical functions in the subsurface are associated with bacteria from new phyla and little studied lineages](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.12930),[**32--Genomic resolution of linkages in carbon, nitrogen, and sulfur cycling among widespread estuary sediment bacteria**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0077-6) and newly reconstructed genomes from current metagenome projects (see Supplementary Table 1 for NCBI accession numbers).

*Accession numbers for all genomes included in this study. Genomes were mined for bacterial-type (P-type) ribosomal proteins L2, L3, L4, L5, L6, L14, L15, L16, L18, L22, L24, S3, S8, S10, S17, and S19 (all single-copy genes), and one representative SSU rRNA gene*

From IMG-M, genomes were sampled such that a single representative for each defined **genus** was selected.
For phyla and candidate phyla lacking full taxonomic definition, every member of the phylum was initially included.
Subsequently, these radiations were sampled to an **approximate genus level** of divergence based on comparison with taxonomically described phyla, thus removing strain- and species-level overlaps.
Finally, initial tree reconstructions identified aberrant long-branch attraction effects placing the Microsporidia, a group of parasitic fungi, with the Korarchaeota.
The Microsporidia are known to contribute long branch attraction artefacts confounding placement of the Eukarya [ref--Covarion Shifts Cause a Long-Branch Attraction Artifact That Unites Microsporidia and Archaebacteria in EF-1α Phylogenies](https://academic.oup.com/mbe/article/21/7/1340/1080413), and were subsequently removed from the analysis.

This study includes 1,011 organisms from **lineages** for which genomes were not previously available.
The organisms were present in samples collected from a shallow aquifer system, a deep subsurface research site in Japan, a salt crust in the Atacama Desert, grassland meadow soil in northern California, a CO2-rich geyser system, and two dolphin mouths.
[Genomes were reconstructed from metagenomes as described previously](https://www.ncbi.nlm.nih.gov/pubmed/26083755).
Genomes were only included if they were estimated to be **>70% complete** based on presence/absence of a suite of 51 single copy genes for Bacteria and 38 single copy genes for Archaea.
Genomes were additionally required to have consistent nucleotide composition and coverage across scaffolds, as determined using the ggkbase binning software (ggkbase.berkeley.edu), and to show consistent placement across both **SSU rRNA and concatenated ribosomal protein phylogenies**.
This contributed marker gene information for 1,011 newly sampled organisms, whose genomes were reconstructed for metabolic analyses to be published separately.

[The **concatenated ribosomal protein alignment** was constructed as described previously](https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-22).
- the 16 ribosomal protein data sets (ribosomal proteins L2, L3, L4, L5, L6, L14, L16, L18, L22, L24, S3, S8, S10, S17 and S19) were aligned independently using [MUSCLE v. 3.8.31](https://academic.oup.com/nar/article/32/5/1792/2380623).
- Alignments were trimmed to remove ambiguously aligned C and N termini as well as columns composed of more than 95% gaps.
- Taxa were removed if their available sequence data represented less than 50% of the expected alignment columns (90% of taxa had more than 80% of the expected alignment columns).
- The 16 alignments were concatenated, forming a final alignment comprising 3,083 genomes and 2,596 amino-acid positions.
- A maximum likelihood tree was constructed using [RAxML v. 8.1.24](https://academic.oup.com/bioinformatics/article/22/21/2688/251208), as implemented on the [CIPRES web server](https://ieeexplore.ieee.org/document/5676129), under the LG plus gamma model of evolution (PROTGAMMALG in the RAxML model section), and with the number of bootstraps automatically determined (MRE-based bootstopping criterion).
- A total of 156 bootstrap replicates were conducted under the rapid bootstrapping algorithm, with 100 sampled to generate proportional support values.
- The full tree inference required **3,840 computational hours** on the **CIPRES** supercomputer.

For a **companion SSU rRNA tree**, an alignment was generated from **all SSU rRNA genes** available from the genomes of the organisms included in the **ribosomal protein data set**.
For organisms with multiple SSU rRNA genes, one representative gene was kept for the analysis, selected randomly.
As **genome sampling** was confined to the **genus level**, we do not anticipate this selection process will have any impact on the resultant tree.
All SSU rRNA genes **longer than 600 bp** were aligned using the SINA alignment algorithm through the SILVA web interface38,39.
The full alignment was stripped of columns containing 95% or more gaps, generating a final alignment containing 1,871 taxa and 1,947 alignment positions.
A maximum likelihood tree was inferred as described for the **concatenated ribosomal protein trees**, with RAxML run using the **GTRCAT model of evolution**.
The **RAxML** inference included the calculation of 300 bootstrap iterations (extended majority rules-based bootstopping criterion), with **100 randomly sampled to determine support values**.

---
other links:
[IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies](https://academic.oup.com/mbe/article/32/1/268/2925592)
[IQ-TREE supports a wide range of evolutionary models](http://www.iqtree.org/)
[Substitution models](http://www.iqtree.org/doc/Substitution-Models)
[Models of DNA evolution](https://en.wikipedia.org/wiki/Models_of_DNA_evolution)
[The RAxML v8.2.X Manual](https://sco.h-its.org/exelixis/resource/download/NewManual.pdf) [2 RAxML](http://evomics.org/learning/phylogenetics/raxml/)
[Maximum Likelihood Inference and Models of DNA Sequence Evolution](http://ib.berkeley.edu/courses/ib200/labs/05/lab05.pdf)
[Selecting models of evolution](https://www.kuleuven.be/aidslab/phylogenybook/firstEdition/Chapter10.pdf)
[Trends in substitution models of molecular evolution](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4620419/)
[Evolutionary models of phylogenetic trees](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691382/pdf/12965036.pdf)
[Introduction of phylogenetic analysis](https://bip.weizmann.ac.il/education/course/introbioinfo/03/lect12/phylogenetics.pdf)
[To include or not to include: The impact of gene filtering on species tree estimation](https://www.biorxiv.org/content/biorxiv/suppl/2017/06/12/149120.DC1/149120-1.pdf)
---

To test the effect of site selection stringency on the inferred phylogenies, we stripped the alignments of columns containing up to 50% gaps (compared with the original trimming of 95% gaps).
For the ribosomal protein alignment, this resulted in a 14% reduction in alignment length (to 2,232 positions) and a 44.6% reduction in computational time (∼2,100 h).
For the SSU rRNA gene alignment, stripping columns with 50% or greater gaps reduced the alignment by 24% (to 1,489 positions) and the computation time by 28%.
In both cases, the topology of the tree with the best likelihood was not changed significantly.
The ribosomal protein resolved a two-domain tree with the Eukarya sibling to the Lokiarcheaota, while the SSU rRNA tree depicts a three-domain tree.
The position of the **CPR (candidate phylum reconstructed)** as deep-branching on the ribosomal protein tree and within the Bacteria on the SSU rRNA tree was also consistent. The alignments and inferred trees under the more stringent gap stripping are available upon request.
