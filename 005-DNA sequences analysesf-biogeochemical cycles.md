Published paper: [Connecting biodiversity and potential functional role in modern euxinic environments by microbial metagenomics](https://www.nature.com/articles/ismej2014254)

# DNA sequences analyses

The functional analyses focused on the three main biogeochemical cycles for this type of lakes, that is, carbon (C), nitrogen (N) and sulfur (S) cycling. The genetic potential of the microbial community was analyzed following the C, N, and S marker genes (KOs) as reported by [Lauro et al. (2011)](https://www.nature.com/articles/ismej2010185) with a few modifications.

We amended this previous rubric by adding the **anaerobic carbon fixation** carried out through the Calvin cycle by Chromatiaceae, and additional genes for polysulfide reduction, nitrate reduction and nitrite oxidation.

In addition, the genes pyruvate:ferredoxin oxidoreductase (porA/B) were not considered as marker genes for fermentation as in [Lauro et al. (2011)](https://www.nature.com/articles/ismej2010185), because they are key genes in the reverse tricarboxylic acid cycle used for carbon fixation by Epsilonproteobacteria abundant in our study lakes ([Campbell and Cary, 2004](https://aem.asm.org/content/70/10/6282); [Takai et al., 2005](https://aem.asm.org/content/71/11/7310)).

Because both sulfide oxidation and dissimilatory sulfate reduction pathways are mediated by the same set of genes (aprA, aprB and dsrA) but are found in different families of bacteria, we assigned metagenomic reads to each pathway according to phylogeny, that is, sulfate reduction for Firmicutes and Deltaproteobacteria reads, and sulfide oxidation for Alphaproteobacteria, Betaproteobacteria, Chlorobiaceae and Chromatiaceae. Finally, for the sulfur-oxidizing Epsilonproteobacteria of the order Campylobacterales we specifically searched for sox genes (coding for thiosulfate oxidation) not currently available in the KEGG database.

Marker genes used in the present work are shown in [Supplementary Table S1](https://media.nature.com/original/nature-assets/ismej/journal/v9/n7/extref/ismej2014254x2.pdf). Hierarchical clustering and heatmap plots were generated with R (R Development Core Team, 2012) using the library ‘seriation’.
