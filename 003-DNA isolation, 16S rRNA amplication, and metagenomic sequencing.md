# DNA isolation, 16S rRNA amplicon, and metagenomic sequencing

Published paper: [A metagenomics roadmap to the uncultured genome diversity in hypersaline soda lake sediments](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0548-7)

## step-1

The **colloidal fraction** of each sediment sample (~ 10% of 50 g) was separated from the course sandy fraction by several short (30–60 s) low-speed (1–2,000 rpm in 50 mL Falcon tubes) centrifugation steps and washed with 1–2 M NaCl solution [detail](http://academics.wellesley.edu/Biology/Concepts/Html/molarsolutions.html). 

## step-2

The **pelleted colloidal sediment fraction** was first subjected to 3 cycles of freezing in liquid nitrogen/thawing, then re-suspended in 0.1 M Tris (pH 8)/10 mM EDTA, and then subjected to harsh bead beating treatment. Next, the samples were incubated with lysozyme (15 mg/mL) for 2 h at 37 °C followed by a SDS (10% w/v) and proteinase K (10 μg/mL) treatment for 30 min. at 45 °C. High molecular weight DNA was isolated using phenol/chloroform extraction, quality-checked, and sequenced as described previously [Metagenomic Insights into the Uncultured Diversity and Physiology of Microbes in Four Hypersaline Soda Lake Brines](https://doi.org/10.3389/fmicb.2016.00211).

## step-3

Direct high-throughput sequencing of the DNA was performed on an Illumina HiSeq 2000 platform to generate 150 b paired-end reads. Amplification of the **V4-V6** region of prokaryote 16S rRNA genes using barcoded **926F-1392R** primers, amplicon purification, quantification, and Roche (454)-sequencing was performed together in a batch with brine samples from the same sampling campaigns. 

## step-4

Barcodes and adapter sequences were removed from de-multiplexed amplicon sequence reads and analyzed with the automated NGS analysis pipeline of the [SILVA](https://academic.oup.com/nar/article/41/D1/D590/1069277) rRNA gene database project (SILVAngs 1.3, database release **version 128**) using default parameters. 
The OTUs (**97% identity**) assigned down to the **genus** level were only considered when they had a **relative abundance ≥ 0.1%** in at least one of the five datasets.




