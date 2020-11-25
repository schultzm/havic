---
title: 'havic: detect Hepatitis A Virus Infection Clusters from clinical cDNA sequences'
tags:
  - python3
  - bioinformatics
  - pipeline
  - ruffus
  - sequence mapping
  - consensus
  - epidemiology
  - genomic
  - molecular
  - amplicon
  - HAVNet
  - HAV
  - outbreaks
  - phylogenomics
authors:
  - name: Mark B. Schultz
    orcid: 0000-0002-7689-6531
    affiliation: 1
  - name: Karolina Mercoulia
    orcid: 0000-0002-3957-666X
    affiliation: 1
  - name: William Pitchers
    orcid: 0000-0003-0385-5939
    affiliation: 1
  - name: Patiyan Andersson
    orcid: 0000-0003-2705-1847
    affiliation: 1
  - name: Marion Easton
    affiliation: 3
  - name: Joy Gregory
    affiliation: 3
  - name: Danita Hennessy
    affiliation: 3
  - name: Linda T. Viberg
    orcid: 0000-0002-4174-1142
    affiliation: 1
  - name: Michelle Sait
    orcid: 0000-0001-8782-6128
    affiliation: 1
  - name: Susan A. Ballard
    orcid: 0000-0002-0096-9474
    affiliation: 1
  - name: Norelle L Sherry
    orcid: 0000-0002-7789-8360
    affiliation: 1
  - name: Lilly Yuen
    orcid: 0000-0003-2934-0167
    affiliation: 2
  - name: Leon Caly
    orcid: 0000-0003-2846-6729
    affiliation: 2
  - name: Mike G. Catton
    orcid: 0000-0002-0297-1260
    affiliation: 2
  - name: Julian D.Druce
    orcid: 0000-0002-3521-3298
    affiliation: 2
  - name: Courtney R. Lane
    affiliation: 1
  - name: Kristy Horan
    orcid: 0000-0003-3960-1519
    affiliation: 1
  - name: Torsten Seemann
    orcid: 0000-0001-6046-610X
    affiliation: 1
  - name: Benjamin P. Howden
    orcid: 0000-0003-0237-1473
    affiliation: 1
  - name: Anders Gonçalves da Silva
    orcid: 0000-0002-2257-8781
    affiliation: 1
  - name: Deborah A. Williamson
    orcid: 0000-0001-7363-6665
    affiliation: 1
affiliations:
  - name: Microbiological Diagnostic Unit – Public Health Laboratory (MDU-PHL)
    index: 1
  - name: Victorian Infectious Diseases Reference Laboratory (VIDRL)
    index: 2
  - name: Department of Health and Human Services (DHHS), Victorian Government, Australia
    index: 3
date: 25 November 2020
bibliography: paper.bib
---

# Introduction

Globally, Hepatitis A Virus (HAV) infects tens of millions of people each year [@RN406].  Though mortality is low, morbidity is high resulting in large economic burden due to frequent hospitalisation of cases during the acute phase of infection.  Transmission of HAV occurs via the faecal-oral pathway, either directly from person-to-person or indirectly through contaminated food and water [@RN398].  HAV has its highest prevalence in low-income countries; however, sporadic outbreaks do occur in high-income countries [@RN406], typically arriving via return travellers and import on frozen foods [e.g., see @RN410; @RN414].  Lifelong immunity to the virus arises after vaccination or infection [@RN416].

Molecular epidemiology using DNA sequencing and comparative genomics to characterise virus outbreaks in real-time is now considered an essential public health measure [@RN407; @RN410; @RN415].  For HAV, whole genome sequencing is not yet the normal practice for outbreak surveillance.  Instead, the gold-standard approach [i.e., the 'HAVNet' protocol] involves the sequencing of a 460 bp cDNA amplicon spanning the VP1/P2A junction.  After sequencing, the amplicon is compared to global databases to make inferences of genotype and to recover putative epidemiological links [@RN316].  As tiled amplicon approaches to viral sequencing are becoming more common [e.g., @RN417], there is and will be a need for software tools that can handle amplicon, whole genome and/or partial genome sequences.  

Diversity of HAV has been well characterised and modern nomenclature assigns HAV to one serotype with six genotypes (I-VI).  Of the six, three genotypes (i.e., I, II and III) infect humans, and three (i.e., IV, V and VI) infect non-human primates.  Phylogenetic analysis conspicuously divides genotypes I, II and III into six subtypes (IA, IB, IIA, IIB, IIIA and IIIB) [@RN416; @RN412].  Pairwise nucleotide divergences between genotypes I to III range from 12.2% to 21.9% (mean 18.2%).  Diversity _within_ genotypes I to III ranges from 0.3% to 6.5% (mean 4.3%). Between subtypes IA and IB mean nucleotide divergence is estimated at 9.3%; between IIA and IIB divergence is approximately 9.6%; and between IIIA and IIIB divergence is around 11.8% [@RN412].  

Though the HAVNet protocol is widely adopted, variations to the method are commonplace [e.g., see @RN407; @RN408].  Public databases (e.g., HAVNet and NCBI) are compiled over many years from myriad laboratories [e.g., see @RN411].  Artefactual nucleotide variations (e.g., low quality, false indels, incompatible orientation) are present.  And HAV consensus sequences do not always co-locate within a single genome target, with databases comprised of sequences from whole genome or partial genome sequences from varying regions [refer to Figure 1 in @RN416].  Distance based pairwise nucleotide comparisons used to infer relatedness of samples are naturally sensitive to these sources of variation.  In this article, we describe our software tool `havic`, which is written to ease the burden of detecting and characterising HAV outbreak clusters during routine epidemiolgical surveillance in the midst of these challenges.

# Statement of utility

`havic` has been developed over a number of years with feedback and feature requests from public health epidemiolgists during routine use of the program for HAV outbreak surveillance.  The software aims to provide repeatable analyses that give actionable and objective results, despite inherent imperfections in HAV sequencing data.  As the HAV genome comprises a single segment, being a positive-sense single-stranded ribonucleic acid (RNA) of only 7.5 kilobases (kb) [@RN375], `havic` runs can be completed relatively quickly on a standard desktop computer.  For larger jobs, the software pipeline has been tested on an HPC system with queing managed by SLURM.  

The pipeline is written in python3, implementing a number of R packages, and pipeline control within `havic` is managed using the Ruffus library.  Run configurations are defined using a yaml-formatted text file.  Installation is performed using `conda`.  

# The workflow

Briefly, outlining the salient steps in our pipeline, `havic`:

- reads DNA sequences in fasta format from one or many files (BioPython)  
- finds and removes duplicates (based on fasta header)  
- maps the sequences to a reference genome (reverse complementing the input sequences as required) using `minimap2` [@RN275]  
- converts the sequence alignment map (SAM) output from `minimap2` to a binary alignment map (BAM) using `samtools` [@RN314] (excluding secondary mappings and unmapped reads)  
- extracts the aligned reads from the BAM as a multiple sequence alignment (MSA) using RSamtools [@RN276]
- trims sequences in the MSA to a user-defined target region (BioPython)
- infers a nucleotide substitution model from the MSA using ModelFinder [@RN323]

Under the best-fit substitution model, `havic`:

- runs IQtree2 [@RN315] to infer the most likely tree using Maximum Likelihood (ML)
- infers support values for branches in the tree using the the ultra-fast bootstrap method [@RN418]

After obtainining the tree, `havic`:

- run `ClusterPicker` [@RN274] to delineate putative infection clusters within the input sample set
- reports infection clusters in textual, tree and graphical format files

# Testing and validation

Example data are pre-packaged with `havic`.  Tests are performed by running `havic test <test_suite>`, where `<test_suite>` is any of `hav_amplicon`, `hav_wgs`.  Two additional test suites – `measles_wgs` and 'hiv_amplicon' – are included for exploration of edge-cases and for future work that will aim to allow iterations of the software to be used for outbreak surveillance in other viruses.  Considering the error rate of Q30 Sanger sequences is 1/1000 bases or 0.001, and based on advice from epidemiolgists working in this domain (and who have been independently able to validate `havic`-detected clusters using contact tracing metadata over a number of years and multi-jurisdictional outbreaks), we have settled on an optimal maximal within infection cluster divergence of 0.01 for our test analyses.  Using the `UFBoot` method in `IQ-Tree2`, we performed our tests using branch supports of greater than or equal to 95%.  The software was developed using a test driven ethos, with python unittests forming the basis of all tests.  As a `havic` run makes inferences from a single best tree inferred under ML, we recommend running the pipeline multiple times to avoid reliance on a sub-optimal tree.  

# Visualisation tools

The results of havic are output to a single folder with the option to summarise the results as images.  Pairwise nucleotide differences may be summarised as a heatmap \autoref{fig:heatmap}

![Figure 1: Pairwise genetic distances and ClusterPicker clusters.\label{fig:heatmap}](https://github.com/schultzm/havic/blob/master/havic/data/_heatmap_SNPs.png?raw=true )

And the alignment may be plotted next to the phylogenetic tree, with the tree tips coloured by infection cluster.  
![Tree](https://github.com/schultzm/havic/blob/master/havic/data/tree_MSA_clusters.png?raw=true "Maximum Likelihood tree with bootstrap support, ClusterPicker clusters, and Multiple Sequence Alignment")

Samples listed under HIGHLIGHT_TIP will be annotated in the final tree plot with a red dot, as shown below.  

![Tree](https://github.com/schultzm/havic/blob/master/havic/data/highlight_tip.png?raw=true "Tip CmvAXJTIqH highlighted as requested under HIGHLIGHT_TIP")


# Acknowledgements

We would like to thank members of the laboratories at Queensland Health Forensic and Scientific Services (QHFSS), DHHS, VIDRL and MDU-PHL.

# References
