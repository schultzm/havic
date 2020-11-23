---
title: 'havic: a software tool for detection of Hepatitis A Virus Infection Clusters from clinical cDNA sequences'
tags:
  - Python3
  - Python
  - bioinformatics
  - pipeline
  - ruffus
  - mapping
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
  - name: Mark B. Schultz^[Corresponding author]
    orcid: 0000-0002-7689-6531
    affiliation: 1
  - name: Karolina Mercoulia
    affiliation: 1
    orcid: 0000-0002-3957-666X
  - name: William Pitchers
    orcid: 0000-0003-0385-5939
    affiliation: 1
  - name: Patiyan Andersson
    orcid: 0000-0003-2705-1847
    affiliation: 1
  - name: Marion Easton
    affiliation: 3
    orcid: NA
  - name: Joy Gregory
    affiliation: 3
    orcid: NA
  - name: Danita Hennessy
    affiliation: 3
    orcid: NA
  - name: Linda T. Viberg
    orcid: 0000-0002-4174-1142
    affiliation: 1
  - name: Michelle Sait
    orcid: 0000-0001-8782-6128
    affiliation: 1
  - name: Susan A. Ballard
    orcid: 0000-0002-0096-9474
    affiliation: 1
  - Norelle L Sherry
    orcid: 0000-0002-7789-8360
    affiliation: 1
  - Lilly Yuen
    orcid: 0000-0003-2934-0167
    affiliation: 2
  - Leon Caly
    orcid: 0000-0003-2846-6729
    affiliation: 2
  - Mike G. Catton
    orcid: 0000-0002-0297-1260
    affiliation: 2
  - Julian D.Druce
    orcid: 0000-0002-3521-3298
    affiliation: 2
  - Courtney R. Lane
    orcid: 
  - name: Kristy Horan
    orcid: 0000-0003-3960-1519
    affiliation: 1
  - name: Torsten Seemann
    orcid: 0000-0001-6046-610X
    affiliation: 1
  - name: Deborah A. Williamson
    orcid: 0000-0001-7363-6665
    affiliation: 1
  - name: Anders Gonçalves da Silva
    orcid: 0000-0002-2257-8781
    affiliation: 1
  - name: Benjamin P. Howden
    orcid: 0000-0003-0237-1473
    affiliation: 1

affiliations:
  - name: Microbiological Diagnostic Unit – Public Health Laboratory
    index: 1
    address: Department of Microbiology and Immunology, University of Melbourne at The Peter Doherty Institute for Infection and Immunity, 792 Elizabeth Street, Melbourne, Victoria, Australia, 3000
  - name: Victorian Infectious Diseases Reference Laboratory
    index: 2
    address: Melbourne Health, The Peter Doherty Institute for Infection and Immunity, 792 Elizabeth Street, Melbourne, Victoria, Australia, 3000
  - name: Department of Health and Human Services, Victorian Government, Australia

date: 21 September 2020
bibliography: paper.bib
csl: harvard-the-university-of-melbourne.csl
---

# Introduction

Globally, Hepatitis A Virus (HAV) infects tens of millions of people each year [@RN406].  Though mortality is low, morbidity is high resulting in large economic burden due to frequent hospitalisation of cases during the acute phase of infection.  Transmission of HAV occurs via the faecal-oral pathway, either directly from person-to-person or indirectly through contaminated food and water [@RN398].  HAV has its highest prevalence in low-income countries; however, sporadic outbreaks do occur in high-income countries [@RN406], typically arriving via return travellers and import on frozen foods [@RN410, @RN414].  Lifelong immunity to the virus arises after vaccination or infection [@RN416].  

Molecular epidemiology using DNA sequencing and comparative genomics is now considered an essential public health measure to characterise virus outbreaks in real-time [@RN407, @RN410, @RN415].  For HAV, whole genome sequencing is not yet the normal practice for outbreak surveillance.  Instead, the gold-standard approach [i.e., the 'HAVNet' protocol @RN316] involves the sequencing of a 460 bp cDNA amplicon spanning the VP1/P2A junction.  The amplicon is then compared to global databases to make inferences of genotype and to recover putative epidemiological links [@RN316].

Using sequence data, genetic diversity of HAV has been well characterised and modern nomenclature assigns HAV to one serotype divided into six genotypes (I-VI).  Three genotypes (I, II, III), comprising six subtypes (IA, IB, IIA, IIB, IIIA and IIIB), are known to infect humans [@RN416].  Genotypes IV, V and VI infect non-human primates.  Genetic divergence between coding regions of non-epidemiologically linked whole genome sequences is up to 12% within HAV.  Pairwise nucleotide divergences between genotypes I to III range from 12.2% to 21.9% (mean, 18.2%).  Diversity within genotypes I to III ranges from 0.3% to 6.5% (mean, 4.3%). Between subtypes IA and IB mean nucleotide divergence is 9.3%; between IIA and IIB divergence is 9.6%; and between IIIA and IIIB divergence is 11.8% [@RN412].  

Though the HAVNet protocol is widely adopted, variations to the method are commonplace [e.g., see @RN408, @RN407].  Public databases (e.g., HAVNet and NCBI) are compiled over many years from myriad laboratories [e.g., see @RN411].  Artefactual nucleotide variations (e.g., low quality, false indels, incompatible orientation) are present.  And HAV consensus sequences do not always co-locate within a single genome target [refer to Figure 1 in @RN416].  Distance based pairwise sequence comparisons used to infer relatedness of samples are naturally sensitive to these sequence variations.  In this article, we describe our software tool `havic`, which is written to detect and characterise HAV outbreaks in the midst of these challenges with the aim of easing the burden of detecting outbreak clusters during routine epidemiological surveillance.  

# Statement of utility

`havic` has been developed over a number of years with feedback and feature requests from public health epidemiolgists during routine use of the program for HAV outbreak surveillance.  The software aims to provide repeatable analyses that give actionable and objective results, despite inherent imperfections in HAV sequencing data.  As the HAV genome comprises a single segment, being a positive-sense single-stranded ribonucleic acid (RNA) of only 7.5 kilobases (kb) [@RN375], `havic` runs can be completed relatively quickly on a standard desktop computer.  For larger jobs, the software pipeline has been tested on a HPC system with queing managed by SLURM.  

The pipeline is written in python3, implementing a number of R packages, with pipeline control within `havic` managed using the Rufus library.  Run configurations are defined using a yaml-formatted text file.  Installation is performed using `conda`.  


# Testing and validation

Example data are pre-packaged with `havic`.  Tests are performed by doing `havic test <test_suite>`, where `<test_suite>` is any of `hav_amplicon`, `hav_wgs`.  Two additional test suites – `measles_wgs` and 'hiv_amplicon' – are included for exploration of edge-cases and for future work in , with the aim of providing a t other viruses.  , .  Previous analyses are difficult to replicate without this innformation (e.g. italy paper, clusterpicker demo).  Q30 error rate is 1/1000 or 0.1 in 100 bases, translating to a distance of 0.001.  

# Figures


# Acknowledgements


# References

