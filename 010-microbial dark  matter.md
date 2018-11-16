pulished paper:

[Insights into the phylogeny and coding potential of microbial dark matter](https://www.nature.com/articles/nature12352) 2013

## many papers cited this one

[SUPPLEMENTARY INFORMATION](https://media.nature.com/original/nature-assets/nature/journal/v499/n7459/extref/nature12352-s1.pdf)

# SAG (single-cell amplification) assembly

The draft genome of all but 13 SAGs was generated at the DOE Joint genome Institute (JGI) using the Illumina technology. An Illumina standard shotgun library was constructed and sequenced using the Illumina HiSeq 2000 platform. All general aspects of library construction and sequencing performed at the JGI can be found at http://www.jgi.doe.gov/.

Raw Illumina sequence data were filtered for known Illumina sequencing and library preparation artifacts and then screened and trimmed according to the k-mers present in the dataset.

High-depth k-mers, presumably derived from MDA amplification bias, cause problems in the assembly, especially if the k-mer depth varies in orders of magnitude for different regions of the genome. Reads representing highly abundant k-mers were removed such that no k-mers with a coverage of more than 30x were present after filtering.

Reads with an average kmer depth of less than 2x were removed.

The following steps were then performed for assembly:
- (1) filtered Illumina reads were assembled using Velvet version 1.1.04 (Zerbino and Birney, 2008). The VelvetOptimiser script (version 2.1.7) was used with default optimization functions (n50 for k-mer choice, total number of base pairs in large contigs for cov_cutoff optimization).
- (2) 1-3 kbp simulated paired end reads were created from Velvet contigs using the wgsim software.
- (3) the normalized Illumina reads were assembled together with simulated read pairs using Allpaths-LG (version 41043) (Gnerre and MacCallum, 2011).
---
Parameters for assembly steps were:
1) VelvetOptimiser (--v --s 51 --e 71 --i 4 --t 1 --o "-ins_length 250 -min_contig_lgth 500")
2) wgsim (-e 0 -1 100 -2 100 -r 0 -R 0 -X 0)
3) Allpaths-LG (prepareAllpathsParams: PHRED_64=1 PLOIDY=1 FRAG_COVERAGE=125 JUMP_COVERAGE=25 LONG_JUMP_COV=50, runAllpathsParams: THREADS=8 RUN=std_pairs TARGETS=standard VAPI_WARN_ONLY=True OVERWRITE=True).
---
*For the remaining 13 SAGs the draft genomes were generated at the JGI using a Roche 454 Genome Sequencer FLX System using Titanium chemistry according to the manufacturer’s protocols (454 Life Sciences, Branford, CT). The hybrid 454/Illumina assemblies steps (1)-(3) were identical.*
---
- (4) Allpaths contigs larger than 1 kbp were shredded into 1 kbp pieces with 200 bp overlaps.
- (5) the Allpaths shreds and raw 454 pyrosequence reads were assembled using the 454 Newbler assembler version 2.4 (Roche). All assemblies are available on the Microbial Dark Matter project website (http://genome.jgi.doe.gov/MDM/)
