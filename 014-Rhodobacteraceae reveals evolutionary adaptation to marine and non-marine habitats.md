published paper: [2017 -- Phylogenomics of Rhodobacteraceae reveals evolutionary adaptation to marine and non-marine habitats](https://www.nature.com/articles/ismej2016198)

# Materials and methods

## 01-Genome-scale phylogenetic analysis

All applied methods are detailed in [**Supplementary File S1**](https://media.nature.com/original/nature-assets/ismej/journal/v11/n6/extref/ismej2016198x2.pdf).

Among the 106 genome-sequenced strains investigated, 13 strains of the Labrenzia/Stappia group taxonomically assigned to Rhodobacteraceae but rather placed within Rhizobiales in 16S rRNA gene analyses [(Pujalte et al., 2014)](https://link.springer.com/referenceworkentry/10.1007%2F978-3-642-30197-1_377) were used as outgroup.

**An extended data set** including 132 genomes was phylogenetically analysed using **Genome BLAST Distance Phylogeny** ([Auch et al., 2006 -- Genome BLAST distance phylogenies inferred from whole plastid and whole mitochondrion genome sequences](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-350); [Meier-Kolthoff et al., 2014 -- Highly parallelized inference of large genome‐based phylogenies](https://onlinelibrary.wiley.com/doi/abs/10.1002/cpe.3112)).

**Digital DNA:DNA hybridization** was used to check all **species affiliations** ([Auch et al., 2010--Standard operating procedure for calculating genome-to-genome distances based on high-scoring segment pairs](http://standardsingenomics.org/content/2/1/142/); [Meier-Kolthoff et al., 2013a--Genome sequence-based species delimitation with confidence intervals and improved distance functions](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-60)).

**Pairwise 16S rRNA gene similarities** ([Meier-Kolthoff et al., 2013b -- When should a DDH experiment be mandatory in microbial taxonomy?](https://link.springer.com/article/10.1007%2Fs00203-013-0888-4)) were determined after extraction with **RNAmmer version 1.2** ([Lagesen and Hallin, 2007](https://academic.oup.com/nar/article/35/9/3100/2401119)).

The **proteome sequences** were phylogenetically investigated using the DSMZ phylogenomics pipeline ([Anderson et al., 2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020237); [Breider et al., 2014](https://www.frontiersin.org/articles/10.3389/fmicb.2014.00416/full); [Frank et al., 2014](http://standardsingenomics.org/content/9/3/914/)).

**Alignments were concatenated** to three main supermatrices:
- (i) ‘core genes’, alignments containing sequences from all proteomes;
- (ii) ‘full’, alignments containing sequences from at least four proteomes;
- (iii) ‘MARE’, the full matrix filtered with that software ([Meusemann et al., 2010--A Phylogenomic Approach to Resolve the Arthropod Tree of Life](https://academic.oup.com/mbe/article/27/11/2451/1113513)).

The **core genes** were further reduced to their 50, 100, 150 and 200 most conserved genes (up to 250 without outgroup).

**Long-branch extraction** ([Siddall and Whiting, 1999](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1096-0031.1999.tb00391.x)) to assess long-branch attraction artefacts ([Bergsten, 2005--A review of long‐branch attraction](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1096-0031.2005.00059.x)) was conducted by removing the outgroup strains, generating the supermatrices anew and rooting the resulting trees with **LSD version 0.2** ([To et al., 2015--Fast Dating Using Least-Squares Criteria and Algorithms](https://academic.oup.com/sysbio/article/65/1/82/2461506)).

ML and maximum parsimony (MP) phylogenetic trees were inferred as described ([Andersson et al., 2011](https://doi.org/10.1111/j.1096-0031.2008.00217.x); [Breider et al., 2014](https://doi.org/10.3389/fmicb.2014.00416); [Frank et al., 2014](https://doi.org/10.4056/sigs.5179110)) but **MP tree** search was conducted with **TNT version 1.1** ([Goloboff et al., 2008](https://doi.org/10.1111/j.1096-0031.2008.00217.x)).

Additionally, best **substitution models** for each gene and ML phylogenies were calculated with **ExaML version 3.0.7** ([Stamatakis and Aberer, 2013](https://ieeexplore.ieee.org/document/6569896)).

Ordinary and partition bootstrapping ([Siddall, 2010](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1096-0031.2009.00295.x)) was conducted with 100 replicates except for full/ML for reasons of running time.

To assess conflict between the genomic data and roseobacter monophyly, site-wise ML and MP scores were calculated from unconstrained and accordingly constrained best trees, optionally summed up per gene, and compared with Wilcoxon and T-tests (R) and, for ML, with the approximately unbiased test ([Shimodaira and Hasegawa, 2001](https://doi.org/10.1093/bioinformatics/17.12.1246)).

Potential conflict with the 16S rRNA gene was measured using constraints derived from the supermatrix trees.

Major sublineages were inferred from the phylogenomic trees non-arbitrarily as the maximally inclusive, maximally and consistently supported subtrees.

## 02-Analysis of character evolution

Phylogenetic correlations between pairs of binary characters ([Pagel, 1994](https://doi.org/10.1098/rspb.1994.0006)) were detected with **BayesTraits version 2.0** ([Pagel et al., 2004](https://doi.org/10.1080/10635150490522232)) in conjunction with the rooted **ML phylogenies**.

**Ratios** of the estimated rates of change were calculated to verify the tendency toward marine or non-marine habitats.

Three distinct genome samplings were used to detect an influence of only partially sequenced genomes; only results stable with respect to topology and genome sampling were considered further.

Evolution of selected genomic characters was visualized using **Mesquite v2.75** ([Maddison and Maddison, 2011](http://www.mesquiteproject.org/)).

Habitat assignments ([Supplementary File S2](https://media.nature.com/original/nature-assets/ismej/journal/v11/n6/extref/ismej2016198x3.xls)) found in the literature only allowed for distinguishing marine and non-marine habitats but this distinction was fully supported by isolation location, physiology and environmental sequencing wherever available ([Supplementary File S3](https://media.nature.com/original/nature-assets/ismej/journal/v11/n6/extref/ismej2016198x4.pdf)).

Habitats with a salt concentration comparable to the sea were considered marine ([Hiraishi and Ueda, 1994](https://doi.org/10.1099/00207713-44-1-15); [Brinkhoff and Muyzer, 1997](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?holding=npg&cmd=Retrieve&db=PubMed&list_uids=9327542&dopt=Abstract)).

The **EnzymeDetector** ([Quester and Schomburg, 2011-- EnzymeDetector: an integrated enzyme function prediction tool and database](https://doi.org/10.1186/1471-2105-12-376)) was used for initial enzyme annotations of the genomes, improved by **strain-specific information** from BRENDA and AMENDA ([**Schomburg et al., 2013**--BRENDA in 2013: integrated reactions, kinetic data, enzyme function data, improved disease classification: New options and contents in BRENDA](https://doi.org/10.1093/nar/gks1049)), BrEPS ([**Bannert et al., 2010**--BrEPS: a flexible and automatic protocol to compute enzyme-specific sequence profiles for functional annotation](https://doi.org/10.1186/1471-2105-11-589)) and BLAST (Altschul et al., 1990) search against UniProt ([**UniProt Consortium, 2013**](https://doi.org/10.1093/nar/gks1068)).

To validate the completeness of each proteome, its proportion of enzymes was calculated (see Supplementary File S2).

Enzymes were mapped on **MetaCyc pathways** (Caspi et al., 2012) as *previously described ([Chang et al., 2015--BRENDA in 2015: exciting developments in its 25th year of existence](https://doi.org/10.1093/nar/gku1068))*.

A pathway with **≥75%** of its enzymes present was initially assumed to be present too; for pathways discussed in detail, this was manually refined considering the enzymes essential for pathway functionality.

The genomic clusters of orthologous groups (COGs) were taken from Integrated Microbial Genomes ([Mavromatis et al., 2009](https://doi.org/10.4056/sigs.632)) Prodigal annotations ([Hyatt et al., 2010](https://doi.org/10.1186/1471-2105-11-119)).

Plasmid replication systems and compatibility groups were determined as described (Petersen et al., 2009, 2011), as were **flagellar gene clusters and flagellar types** ([Frank et al., 2015a--Ocean’s twelve: flagellar and biofilm chromids in the multipartite genome of Marinovum algicola DG898 exemplify functional compartmentalization in Proteobacteria](https://doi.org/10.1111/1462-2920.12947)).

For details, as well as for tests for phylogenetic inertsia (Diniz-Filho et al., 1998) and quantification of oligotrophy (Lauro et al., 2009), see Supplementary File S1.
