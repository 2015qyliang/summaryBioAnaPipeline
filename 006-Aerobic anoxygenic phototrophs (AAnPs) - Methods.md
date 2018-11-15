Published paper: [Potential for primary productivity in a globally-distributed bacterial phototroph](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018677/#CR12)

# Aerobic anoxygenic phototrophs (AAnPs) - Methods

A collection of 3,655 non-redundant MAGs generated from several studies using the Tara Oceans metagenomic dataset ([Tully, Graham, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5769542/), n = 2,631; [Tully, Sachdeva, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5507172/), n = 201; [Delmont et al., 2017](https://www.biorxiv.org/content/early/2017/04/23/129791), n = 688) and from the Red Sea ([Haroon et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4932879/), n = 135) had putative DNA coding sequences (CDSs) predicted using Prodigal [(Hyatt et al., 2012)](https://academic.oup.com/bioinformatics/article/28/17/2223/246063) [GitHub/Prodigal](https://github.com/hyattpd/Prodigal) (v.2.6.2; -m -p meta) and annotated by the KEGG database [(Kanehisa, Sato, Kawashima, et al., 2016)](https://academic.oup.com/nar/article/44/D1/D457/2502600) using BlastKOALA [(Kanehisa, Sato and Morishima, 2016)](https://doi.org/10.1016/j.jmb.2015.11.006) (taxonomy group, Prokaryotes; database, genus_prokaryotes + family_eukaryotes). Assessment of pathways and metabolisms of interest were determined using the script [KEGG-decoder.py](www.github.com/bjtully/BioData/tree/master/KEGGDecoder).

Genomes were screened for the predicted presence of genes assigned as the M and L subunits of **type-II photochemical reaction center** (PufML) and the large and small units of **ribulose-1,5-bisphosphate carboxylase** (RbcLS).

After identification of the AAnP genomes of interest, genomes were subjected to manual assessment of quality.

- For genomes from [Tully, Graham, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5769542/) [Tully, Sachdeva, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5507172/), read coverage and DNA compositional data was utilized to bin additional putative contigs into the genomes of interest using CONCOCT [(Alneberg et al., 2014)](http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/CONCOCT.html) [GitHub](https://github.com/BinPro/CONCOCT) (v.0.4.1; parameters: -c 800 -I 500) on contigs >5kb from each province with an AAnP genome. To improve completion estximates, overlapping [CONCOCT](https://github.com/BinPro/CONCOCT) and [BinSanity](https://github.com/edgraham/BinSanity) [(Graham et al., 2017)](https://peerj.com/articles/3035/) bins were visualized using [Anvi’o (Eren et al., 2015)](https://peerj.com/articles/1319/) (v.2.1.0) and manually refined to improve genome completion and minimize contamination estimates.

- Genomes from [Delmont et al., 2017](https://www.biorxiv.org/content/early/2017/04/23/129791) were visualized in Anvi’o and manually curated based on DNA composition (%G+C and tetranucleotide frequencies).

To determine phylogeny, a set of **31 single copy marker genes** [(Santos and Ochman, 2004)](https://doi.org/10.1111/j.1462-2920.2004.00617.x) were identified in the MAGs of interest that consisted of **ribosomal proteins**, often used for phylogenomic analysis [(Hug et al., 2016)](https://www.nature.com/articles/nmicrobiol201648), and proteins essential for cellular metabolism [(Supplemental Information 2)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018677/bin/41396_2018_91_MOESM5_ESM.xlsx).

A total of 1,099 reference genomes [(Supplemental Information 3)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018677/bin/41396_2018_91_MOESM6_ESM.xlsx) accessed from NCBI [GenBank (Benson et al., 2005)](https://academic.oup.com/nar/article/35/suppl_1/D21/1116856) using [HMMER (Finn et al., 2011)](https://academic.oup.com/nar/article/39/suppl_2/W29/2506513) (v3.1b2; hmmsearch -E 1e-5) and hidden Markov models collected from the Pfam (Bateman et al., 2002) and TIGRfam (Haft et al., 2003) databases.

Genomes with ≥17 markers were used for phylogenetic placement. Each individual marker gene was aligned using MUSCLE (Edgar, 2004) (v3.8.31; parameter: -maxiters 32), trimmed using [TrimAL (Capella-Gutiérrez et al., 2009)](https://academic.oup.com/bioinformatics/article/25/15/1972/213148) (v.1.2rev59; parameter: -automated1), manually assessed, and concatenated.

A maximum likelihood tree was generated using [FastTree (Price et al., 2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/) (v.2.1.10; parameters: -lg -gamma).

None of the MAGs possessed a full-length 16S rRNA gene sequence. A 90 bp 16S rRNA gene fragment detected in NAT185 was compared to sequences in the NCBI NT database (Benson et al., 2005) using BLAST with default parameters (Altschul et al., 1997).

Nine full-length 16S rRNA gene sequences were identified that had 100% nucleotide identity over the full length of the NAT185 gene fragment. As the taxonomic information that can be confidently achieved from such a small gene fragment is tenuous, the full-length sequences with 100% match were aligned using the SINA web portal aligner (Pruesse et al., 2012) and added to the non-redundant 16S rRNA gene database (SSURef128 NR) in ARB (Ludwig et al., 2004; v6.0.3) using the Parsimony Quick tool (default parameters).

These nine sequences clustered with a clade of uncultured organisms within the Alphaproteobacteria associated with the Rhodobacteraceae [(Supplemental Information 4)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018677/bin/41396_2018_91_MOESM7_ESM.txt).

Sequences for **PufM** were collected from previously described lineages (Yutin et al., 2005; 2007) and bacterial artificial chromosomes (Béjá et al., 2002) and from genomes in Integrated Microbial Genomes (Markowitz et al., 2006) based on their KEGG Ontology annotation (K08929; Supplemental Information 5 and 6).

Global Ocean Survey (GOS) assemblies (Venter et al., 2004) with predicted CDS (as above) were searched using the collected PufM sequences with a DIAMOND BLASTP search (Buchfink et al., 2014) (v.0.8.36.98; default settings; Supplemental Information 2 and 3).

RbcL sequences were collected from previously described lineages (Tabita et al., 2007; Supplemental Information 7 and 8).

Two separate phylogenetic trees were constructed (as above).
The relative contribution of the ‘Ca. Luxescamonaceae’ genomes to the overall measured planktonic communities was calculated by two different methods in their respective manuscripts.

Briefly, for MAGs described in Tully et al. (2017), an approximate relative abundance was determined based on the total the bacterial and archaeal community signature as detected by a set of 100 single-copy gene markers (Albertsen et al. 2013). This was achieved by searching for all copies of the gene markers in the total, assembled metagenomes from Tully et al. (2017) and identifying which of these markers could be assigned to a MAG.

A length normalized relative abundance value was calculated for each MAG in each sample (Tully, Graham, et al., 2017). Single-copy marker genes for the Delmont et al. (2017) MAGs were detected in the Tully et al. (2017) assembled metagenomes from the Mediterranean and North Atlantic using BLAST (Altschul et al., 1997) and used to calculate relative abundance.

Briefly, for MAGs described in Delmont et al. (2017), reads were recruited to the contigs of each MAG and used to calculate the relative fraction of aligned reads from the total metagenome (Delmont et al., 2017).
