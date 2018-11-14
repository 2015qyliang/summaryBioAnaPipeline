Published paper: [A metagenomics roadmap to the uncultured genome diversity in hypersaline soda lake sediments](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0548-7)

# Processing metagenomics reads, assembly, binning, post-binning and Detailed genome analyses

## step-1
- Metagenomic raw reads were **quality trimmed** using [**Sickle**](https://github.com/najoshi/sickle) (version 1.33), and only reads ≥ 21 b were retained.
- The prokaryotic community structure at taxonomic top levels was extrapolated from ten million randomly sampled singletons from each dataset. Candidate **16S** rRNA fragments > 90 b were identified [63???](https://www.nature.com/articles/srep00135) and compared against the **SILVA SSU database** 128 (<u>==**blastn**</u>==, min. length 90, min. identity 80%, e value 1e-5). To verify that the microbial community composition was indeed mostly prokaryotic, we did a more general screening of the metagenomics reads that identified also candidate **18S** rRNA fragments > 90 b (see Additional file 1: Tables S4-S5).
- The complete trimmed read sets were assembled into contigs ≥ 1 kb with [**MEGAHIT**](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) (https://github.com/voutcn/megahit) using paired-end mode, k min = 21, k max = 131, k step = 10. 
- Genes were predicted using [Prodigal](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) and RNAs with [rna_hmm3](https://doi.org/10.1002/9781118010518.ch44) (http://weizhong-lab.ucsd.edu/meta_rna/) and [tRNAscan-SE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC146525/). Assembled 16S rRNA sequences were compared to a manually curated version from the SILVA SSU database (e value ≥ 1e-5). 
**rRNA finder tools** (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4646959/)(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5372776/)
- Predicted protein sequences were annotated against KEGG with [**GhostKOALA**](https://doi.org/10.1016/j.jmb.2015.11.006) (genus_prokaryotes + family_eukaryotes + viruses) [68]. 
- Marker genes for central metabolic pathways and key environmental element transformations were identified based on K number assignments  [15-Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system](https://www.nature.com/articles/ncomms13219) [69-An integrative study of a meromictic lake ecosystem in Antarctica](https://www.nature.com/articles/ismej2010185) [70-Potential for microbial H2 and metal transformations associated with novel bacteria and archaea in deep terrestrial subsurface sediments](https://www.nature.com/articles/ismej201739) [71-Connecting biodiversity and potential functional role in modern euxinic environments by microbial metagenomics](https://www.nature.com/articles/ismej2014254)

## step-2

Contigs ≥ 2.5 kb were binned with [**METABAT** ](https://peerj.com/articles/1165/) (superspecific mode) based on differential coverage values obtained by mapping all five trimmed readsets to all five contig sets with **==[Bowtie2 ](https://www.nature.com/articles/nmeth.1923)==**. 

The bins were subjected to post-binning ![S13](C:\Users\asus\Desktop\S13.png). 
Bins were assessed with [**==lineage-specific single copy genes==**](http://www.cnki.net/kcms/doi/10.16288/j.yczz.14-392.html) using [CheckM](https://genome.cshlp.org/content/25/7/1043) and further processed with the metagenomics workflow in [**Anvi’o**](https://peerj.com/articles/1319/#p-9). 

Since Candidate Phyla Radiant (CPR) is not included in the **CheckM** reference trees and are likely to have low-genome completeness, we used an existing training file of 797 CPR genomes to identify putative CPR bins [**Important**](http://merenlab.org/2016/04/17/predicting-CPR-Genomes/). Bins with CheckM-completeness ≥ 50% (884 out of 1778) and an additional four CPR bins were further processed. 

Coding sequences were annotated for ****<u>==**taxonomy against NCBI-nr==</u>** (July, 2017) with [**USEARCH**](https://academic.oup.com/bioinformatics/article/26/19/2460/230188) to verify that most hits in each bin were to prokaryotic references. 

Phage or viral contigs were manually removed. 

Genome contamination (redundancy) was estimated based on marker sets of universal single copy genes identified for  [**Bacteria**](http://www.pnas.org/content/110/14/5540) and  [**Archaea**](https://www.nature.com/articles/nature12352)  **==????==**  as implemented in **==Anvi’o==** 

Genome coverage was obtained by mapping trimmed reads with [**BBMap** ](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) v36.x (kfilter 31, subfilter 15, maxindel 80). Bins with ≥ 5% redundancy were further refined with **==Anvi’o==** using circle phylograms (guide trees tnf-cov: euclidian ward) and scanned again for CPR. Post-binning resulted in a total of 2499 **<u>metagenome-assembled genomes (MAGs)</u>** of which 871 were either medium-quality genome drafts (**CheckM** estimated completeness ≥ 50% and contamination ≤ 10% [Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea](https://www.nature.com/articles/nbt.3893) or lower quality draft genomes from CPR.

## step-3

Phylogeny of the MAGs was assessed based on 16 **single-copy ==ribosomal proteins==** and representative reference genomes of major prokaryote lineages across the tree of life [A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648). 

Individual ribosomal proteins in our MAGs were identified by **==K number assignments==**. Only ribosomal proteins ==≥ 80 aa== were considered. 

Initial maximum-likelihood (ML) trees were constructed to determine which organisms belonged to the Archaea, Bacteria, or CPR with [**FastTree** 2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) (WAG + CAT). Final separate trees for the three distant evolutionary groups were constructed in the same manner. 

Each **ribosomal protein set** was aligned separately with [**MAFFT**](https://academic.oup.com/mbe/article/30/4/772/1073398)  (v7.055b, − auto) and concatenated only if a MAG encoded at least 8 out of 16 proteins. For all trees, a 100× posterior bootstraps analysis was performed. 

Phylogenetic trees were visualized together with genome statistics and abundance information using [iTOL](https://academic.oup.com/nar/article/44/W1/W242/2499315). 

We **==cross-checked==** the taxonomic assignments based on the phylogeny of the ribosomal protein cassette with the top hit contig annotations against NCBI-nr and with the reference lineage obtained with CheckM. Lastly, we manually corrected the MAGs for misplaced 16S rRNA genes. 

The final trees presented in the manuscript were redrawn using [FigTree v1.4.3](http://tree.bio.ed.ac.uk/software/figtree/). 

# Detailed genome analyses

CPR MAGs were re-annotated more thoroughly: genes were predicted with [Prokka](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)
- functional predictions were performed by running [InterProScan 5 ](https://academic.oup.com/bioinformatics/article/30/9/1236/237988) **locally** on the supplied COG, CDD, TIGRFAMs, HAMAP, Pfam, and SMART databases. 
- [BLAST Koala](https://doi.org/10.1016/j.jmb.2015.11.006) was used for KEGG pathway predictions. 
- To find putative **carbohydrate-active enzymes** in all final MAGs, we used the web-resource [dbCAN](https://academic.oup.com/nar/article/40/W1/W445/1076481) to annotate all predicted proteins ≥ 80 aa against [CAZy](https://academic.oup.com/nar/article/37/suppl_1/D233/1003505).

To **identify the ==top ten== abundant MAGs** from each respective dataset, ten million randomly sampled singletons were mapped onto each MAG with a ==cut-off of 95% identity== in minimum of 50 bases. Coverage values were additionally **normalized** for genome size and expressed as reads per kilobase of sequence per gigabase of mapped reads [RPKG](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0611-7). 
A positive score (from 871 to 1) was assigned to each MAG according to the **ranking** of the summed RPKG of MAGs in the high-salinity datasets (B1Sed10 and T1Sed) and a negative score according to the ranking of the summed RPKGs in the moderate salinity datasets (CSSed10, CSSed11, T3Sed10). Both scores were summed to get a **“==salinity preference score==”** with MAGs recruiting preferably from <u>high salinity datasets on the positive end, moderate salinity datasets in the negative end,</u> and those without preference in the middle.

We determined **species delineation** for the most abundant MAGs and their closest reference genomes (NCBI-nr) by Average Nucleotide Identity (ANI) and [conserved DNA-matrices](http://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.64483-0), as follows: ANI ≥ 95%, conDNA ≥ 69% = same species, ANI ≥ 95%, condDNA < 69% = might be same species, ANI < 95%, condDNA < 69% = different species. 
Single gene trees based on maximum likelihood were constructed with untrimmed alignments (MAFFT, L-INS-i model) and FastTree 2 (WAG + CAT, increased accuracy, -spr4 -mlacc 2 -slownni) using 100× bootstraps. 
References were pulled from [eggNOG (v4.5.1)](https://academic.oup.com/nar/article/44/D1/D286/2503059) and supplemented with sequences from NCBI-nr or refined according to [7](https://www.frontiersin.org/articles/10.3389/fmicb.2016.00211/full), [33](https://www.nature.com/articles/ismej201653), [46](http://www.pnas.org/content/115/6/E1166), [92](http://www.pnas.org/content/115/6/E1166), [93](https://febs.onlinelibrary.wiley.com/doi/full/10.1111/j.1742-4658.2012.08811.x), [94](https://mmbr.asm.org/content/71/4/576). 
The curated MAGs were scanned for the presence of **rhodopsin sequences** [1](https://en.wikipedia.org/wiki/Microbial_rhodopsin) [2](https://mmbr.asm.org/content/80/4/929) [3](http://sci-hub.tw/10.1093/database/bav080) with the [hmmsearch software](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195) and a profile hidden Markov model (HMM) of the bacteriorhodopsin-like protein family (Pfam accession number PF01036). The identified sequences with significant similarity were aligned together with a curated database composed of a collection of type-1 rhodopsins, using [MAFFT (L-INS-i accuracy model)](https://academic.oup.com/mbe/article/30/4/772/1073398). This protein alignment was further utilized to construct a maximum likelihood tree with 100× bootstrap with [FastTree 2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) . 
All other genes were identified using the KEGG annotation

- other links
[Bacteriorhodopsin](https://en.wikipedia.org/wiki/Bacteriorhodopsin)
[RuBisCO / Ribulose-1,5-bisphosphate carboxylase/oxygenase](https://en.wikipedia.org/wiki/RuBisCO)
