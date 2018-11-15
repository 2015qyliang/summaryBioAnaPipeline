published paper: [Diverse Marinimicrobia bacteria may mediate coupled biogeochemical cycles along eco-thermodynamic gradients](https://www.nature.com/articles/s41467-017-01376-9)
---
**Marinimicrobia** es un filo candidato de bacterias recientemente propuesto, previamente conocido como **SAR406**, MGA o Marine Group A. Se las ha encontrado principalmente a grandes profundidades tales como el abismo de Challenger, la fosa de las Marianas y la fosa de Puerto Rico. Este filo tiene una baja representación en muestras pelágicas no profundas y alta abundancia en las muestras profundas. Aunque que a menudo a estas bacterias se las relaciona con ambientes de bajo nivel de oxígeno disuelto, se conoce poco sobre su ecología y funciones metabólicas.1​ Marinimicrobia forma parte del grupo FCB junto a otros filos de bacterias relacionados.2​

r1:[Identification of Free-Living and Particle-Associated Microbial Communities Present in Hadal Regions of the Mariana Trench](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4860528/)
r2:[Carbon and Sulfur Cycling below the Chemocline in a Meromictic Lake and the Identification of a Novel Taxonomic Lineage in the FCB Superphylum, Candidatus Aegiribacteria](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4846661/)
---
---
#  Single-cell amplified genomes (SAGs) collection, sequencing, assembly, and decontamination

- Replicate 1-ml aliquots of sea water collected for single-cell analyses were cryopreserved with [**6% glycine betaine (Sigma-Aldrich)**](https://www.sigmaaldrich.com/catalog/product/sigma/b7045?lang=en&region=SG) *why use it??*, frozen on dry ice and stored at −80 °C.
- Single-cell sorting, whole-genome amplification, real-time PCR screens, and PCR product sequence analyses were performed at [the Bigelow Laboratory for Ocean Sciences Single Cell Genomics Center](www.bigelow.org/scgc), as described by [Stepanauskas and Sieracki--Matching phylogeny and metabolism in the uncultured marine bacteria, one cell at a time](http://www.pnas.org/content/104/21/9052)
- SAGs from the NESAP were generated at the DOE Joint Genome Institute (JGI) using the Illumina platform as described in [Rinke et al--Insights into the phylogeny and coding potential of microbial dark matter](https://www.nature.com/articles/nature12352)
- SAGs from Saanich Inlet were sequenced at the Genome Sciences Centre, Vancouver BC, Canada, as described in [Roux et al--Ecology and evolution of viruses infecting uncultivated SUP05 bacteria as revealed by single-cell- and meta-genomics](https://elifesciences.org/articles/03125)
- All SAGs were assembled at JGI as described

The following steps were performed for SAG assembly:
- 1) filtered Illumina reads were assembled using [**Velvet** version 1.1.0437](https://genome.cshlp.org/content/18/5/821) using the VelvetOptimiser script (version 2.1.7) with parameters: (--v --s 51 --e 71 --i 4 --t 1 --o “-ins_length 250 -min_contig_lgth 500”)
- 2) **wgsim** (-e 0 −1 100 −2 100 -r 0 -R 0 -X 0)
- 3) **Allpaths-LG** (prepareAllpathsParams: PHRED_64 = 1 PLOIDY = 1 FRAG_COVERAGE = 125 JUMP_COVERAGE = 25 LONG_JUMP_COV = 50, runAllpathsParams: THREADS = 8 RUN = std_pairs TARGETS = standard VAPI_WARN_ONLY = True OVERWRITE = True)

SAG **prediction analysis and functional annotation** was performed within the [Integrated Microbial Genomes (IMG) platform](https://academic.oup.com/nar/article/40/D1/D115/2902777) (http://img.jgi.doe.gov) developed by the Joint Genome Institute, Walnut Creek, CA, USA

# Construction and validation of population genome bins

Marinimicrobia population genome bins were constructed by identifying metagenomic contigs from Saanich Inlet, and NESAP metagenomes mapping to specific SAG(s) using a supervised binning method based in part on methodologies developed by [Dodsworth et al.-- Single-cell and metagenomic analyses indicate a fermentative and saccharolytic lifestyle for members of the OP9 lineage](https://www.nature.com/articles/ncomms2884) in the construction of OP9 population genome bins.
Initially, determination of membership of individual SAGs to SAG-clusters making up a given phylogenetic clade was conducted.
SAG [**tetranucleotide frequencies**](https://sci-hub.tw/10.1002/elps.1150190412) (**this point  may be a big thing**) were then calculated and converted to z-scores with TETRA (http://www.megx.net/tetra) [44](https://www.ncbi.nlm.nih.gov/pubmed/15305919?dopt=Abstract),[45](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-163).  [Tetra-Nucleotide Analysis (TNA)](http://help.ezbiocloud.net/tetra-nucleotide-analysis-tna/)
[**Z-scores**](https://en.wikipedia.org/wiki/Standard_score) were reduced to three dimensions with *principal component analysis (PCA)* using [PRIMER v6.1.1346](http://updates.primer-e.com/primer7/manuals/Getting_started_with_PRIMER_7.pdf) and *hierarchical cluster analysis* of the z-score PCA with [**Euclidian distance**](https://en.wikipedia.org/wiki/Euclidean_distance) (also performed in PRIMER) was carried out to generate SAG-clusters.
These SAG-clusters reflected phylogenetic placement of the SAGs by SSU rRNA gene analysis.
For construction of population genome bins, metagenomic contigs from NESAP and SI data sets were aligned to SAG contigs with>95% nucleotide identity using **BLAST** and a minimum of 5 kilobase pairs alignment length, Tetranucleotide frequencies of all metagenomic contigs passing this identity and length threshold were calculated and converted to z-scores.
SAG-supervised binning as described in [Dodsworth et al.](https://www.nature.com/articles/ncomms2884) using linear discriminant analysis was carried out using all z-scores with the SAG-bins as training data to classify the metagenomic conigs as making up a given population-genome bin.

Individual SAGs and population genome bins were analyzed for completeness and strain heterogeneity using [CheckM v1.0.554](https://genome.cshlp.org/content/25/7/1043).
Specifically, the **lineage_wf workflow** was used with default parameters.
The **lineage_wf workflow** includes determination of the probable phylogenetic lineage based on detected marker genes.
The determined lineage then dictates the sets of marker genes that is most relevant for estimating a given genome’s completeness and other statistics.
The strain heterogeneity metric is highly informative for population genome bins as it is essentially the **average amino-acid identity** for pairwise comparisons of the (lineage appropriate) redundant single-copy marker genes within a population genome bin (Supplementary Data 5).
For population genome bins the higher the strain heterogeneity value, the more similar the amino acid identity of the redundant maker genes indicating the sequences in the bin originate from a closely related, if not identical, phylogenetic source.

# Marinimicrobia genome **streamlining**

Gene-coding bases and COG-based gene redundancy shown in [Supplementary Fig. 1A, B](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-017-01376-9/MediaObjects/41467_2017_1376_MOESM1_ESM.pdf) were calculated using cluster of orthologous group (COG)-based genome redundancy as described in [Rinke et al.17](https://www.nature.com/articles/nature12352). Each gene’s COG category was predicted through the JGI IMG pipeline. COG redundancy was calculated by averaging the occurrence of each COG in the genome. The percentage of gene-coding bases was calculated by dividing the number of bases contributing to protein- and RNA-coding genes by the total genome size.
For SAGs, the length of the assembled genome was used rather than the estimated genome size.

# Annotation and identification of metabolic genes of interest

Genes of interest were identified in the SAGs and in IMG/M (https://img.jgi.doe.gov/cgi-bin/m/main.cgi) for the metagenomic contigs which made up the population genome bins.
Contigs making up Marinimicrobia population genome bins were run through [**MetaPathways 2.5**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-202) [Metabolic pathways for the whole community](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-619) to annotate **open reading frames (ORFs)** and **reconstruct metabolic pathways**.
As the population genome bins were constructed from multiple metagenomes they contained redundant sequence information, **BLASTp (amino-acid identity cutoff >75%)** was used to identify all copies of a given gene of interest in each population genome bin, which was then used in gene model validation and expression mapping.

# Gene expression mapping

*Metatranscriptomes* from three time points in [Saanich Inlet time series -- A compendium of multi-omic sequence information from the Saanich Inlet water column](https://www.nature.com/articles/sdata2017160) were used to investigate changes in gene expression along water column redox gradients over time for selected ORFs involved in energy metabolism and electron shuttling.
Quality controlled reads from metatranscriptomes were mapped to identified ORFs of interest using [bwa –mem](https://academic.oup.com/bioinformatics/article/25/14/1754/225615) and **reads per kilobase per million mapped (RPKM)** per ORF was calculated using RPKM calculation in [**MetaPathways 2.5**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-202).
For each population genome bin RPKM values for a given sample were summed for ORFs with the same functional annotation to yield an RPKM for a given functional gene.
For other taxonomic groups in Saanich Inlet shown *in Supplementary Fig. 6B*, genes were identified by sequence alignment searches of Saanich Inlet metatranscriptomes (bioSample indicated above) assembled and conceptually translated using **BLASTp** against selected **nitrogen and sulfur cycling genes** from [**Hawley et al**](http://www.pnas.org/content/111/31/11395) and RPKM values calculated as described above.

# **Global distribution** and expression of *nosZ*

Further analysis was carried out to determine the global distribution of Marinimicrobia nosZ in 594 metagenomes.
The nosZ nucleotide sequences from SHBH1141 and ZA3312c, which exhibited a **65% nucleotide identity** to each other by *BLAST*, were clustered at **95% identity** using the *USEARCH cluster* fast algorithm, resulting in three clusters, two SHBH1141 and one ZA3312c.
Nucleotide sequence alignment was carried out using [FAST](https://ieeexplore.ieee.org/document/7758120) [github](https://github.com/hallamlab/FAST), with parameters of **>80%** nucleotide identity and **>60 bp alignment length** against 594 metagenomes.
- For Saanich Inlet and NESAP data sets, abundance of nosZ in a given metagenome or metatranscriptome was determined by **summing the RPKM value for ORF hits** to either SHBH1141 or ZA3312c for a given metagenome or metatranscriptome.
- For 454 sequenced15,42 metagenomes and metatranscriptomes (Peru42 and ETSP15), the number of reads which hit to either SHBH1141 or ZA3312c were summed for a given metagenome.
- For the TARA Oceans data set, the number of genes identified in an assembled metagenome was summed.
- Metatranscriptomic data for Tara was unavailable at this time.
