published paper: [**2015 ISME** - Connecting biodiversity and potential functional role in modern euxinic environments by microbial metagenomics](https://www.nature.com/articles/ismej2014254)

# DNA sequences analyses

## 01-Sequencing  & annotation

Identical reads were removed using [CD-HIT (Li and Godzik, 2006)](https://academic.oup.com/bioinformatics/article/22/13/1658/194225).

Annotation of metagenomic reads was conducted through the [**JCVI prokaryotic annotation pipeline** (Tanenbaum et al., 2010)](http://standardsingenomics.org/content/2/2/229/) using Uniref100, PFAM, TIGRfam and KEGG (Kyoto Encyclopedia of Genes and Genomes) Orthologs (KO) databases for** taxonomic and functional annotation**.

[JCVI Metagenomics reports](http://jcvi.org/metarep) were used for analysis and comparative metagenomics [(Goll et al., 2010)](https://academic.oup.com/bioinformatics/article/26/20/2631/193509).

**KO annotation** was used for **functional analysis** and **KO counts** were **normalized** according to the length of the read and the length of the target gene [(Sharon et al., 2009)](https://link.springer.com/chapter/10.1007/978-3-642-02008-7_35) [presentation](https://pdfs.semanticscholar.org/presentation/50b5/a0d06bdc849b5f324c2326e6d0d9fe001345.pdf).

The communities and functional profiles found in each size fraction were highly similar (Supplementary Figure S1) and, therefore, we **pooled all reads after normalizing** for sequencing depth for subsequent analyses, which allows for a better comparison of metagenomes.

## 02-Functional analyses for carbon (C), nitrogen (N) and sulfur (S) cycling

The functional analyses focused on the three main biogeochemical cycles for this type of lakes, that is, **carbon (C), nitrogen (N) and sulfur (S) cycling**.

The genetic potential of the microbial community was analyzed following the **C, N, and S marker genes (KOs)** as reported by [Lauro et al. (2011)](https://www.nature.com/articles/ismej2010185) with a few modifications.

We amended this previous rubric by *adding* the **anaerobic carbon fixation** carried out through the Calvin cycle by Chromatiaceae, and additional genes for **polysulfide reduction, nitrate reduction and nitrite oxidation**.

In addition, the genes **pyruvate:ferredoxin oxidoreductase (porA/B)** were not considered as marker genes for fermentation as in [Lauro et al. (2011)](https://www.nature.com/articles/ismej2010185), because they are key genes in the **reverse tricarboxylic acid cycle** used for **carbon fixation** by Epsilonproteobacteria abundant in our study lakes ([Campbell and Cary, 2004](https://aem.asm.org/content/70/10/6282); [Takai et al., 2005](https://aem.asm.org/content/71/11/7310)).

Because both **sulfide oxidation** and **dissimilatory sulfate reduction pathways** are mediated by the same set of genes *(aprA, aprB and dsrA)* but are found in different families of bacteria, we **assigned metagenomic reads to each pathway according to phylogeny**, that is,
- sulfate reduction for Firmicutes and Deltaproteobacteria reads, and
- sulfide oxidation for Alphaproteobacteria, Betaproteobacteria, Chlorobiaceae and Chromatiaceae.

Finally, for the **sulfur-oxidizing** Epsilonproteobacteria of the order Campylobacterales we specifically searched for sox genes (coding for thiosulfate oxidation) not currently available in the KEGG database.

**Marker genes used in the present work are shown in** [Supplementary Table S1](https://media.nature.com/original/nature-assets/ismej/journal/v9/n7/extref/ismej2014254x2.pdf).

Hierarchical clustering and heatmap plots were generated with R using the library **‘seriation’**.

Metagenomic data have been *deposited* at [**CAMERA** (Sun et al., 2011)](https://academic.oup.com/nar/article/39/suppl_1/D546/2506708) under accession number CAM_P_0001174.
