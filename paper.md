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
  - molecular epidemiology
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
  - name: Danita Hennessy


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

Globally, Hepatitis A Virus (HAV) infects tens of millions of people each year [@RN406].  Though mortality is low morbidity is high resulting in large economic burden due to frequent hospitalisation of cases during the acute phase of infection.  HAV is transmitted via the faecal-oral pathway directly from person-to-person or through contaminated food and water [@RN398].  The virus has its highest prevalence in low-income countries; however, sporadic outbreaks do occur in high-income countries [@RN406] typically arriving via return travellers and import on frozen foods [@RN410].  

Molecular epidemiology using DNA sequencing and comparative genomics is now considered an essential public health measure to characterise virus outbreaks in real-time [@RN407, @RN410, ].  For HAV, the gold standard genotyping protocol involves sequencing a 460 bp cDNA amplicon of the VP1/P2A junction and comparing this sequence to a global database [see HAVNET @RN316 and NCBI].  Using sequence data, genetic diversity of HAV has been well characterised with one known serotype and seven genotypes (I-VII).  Four of these genotypes (I, II, III and VII), comprising six subtypes (IA, IB, II, IIIA, IIIB, VII), are known to infect humans [@RN404, @RN403].  

Though the HAVNET protocol is widely adopted, variations to the protocol are commonplace [e.g., see @RN408 and @RN407].  Public databases are compiled over many years from myriad laboratories [e.g., see @RN411].  Artefactual nucleotide variation (e.g., low quality, false indels, incompatible orientation) is present.  As a result, HAV consensus sequences do not always co-locate within a single genome target.  Distance based pairwise sequence comparisons used to infer relatedness of samples is thus sensitive to these variations but comparison of large numbers of sequences is often required to resolve outbreaks.  In this article, we describe our software tool `havic`, which is written to characterise HAV outbreaks in the midst of these challenges.  

# Statement of utility

`havic` has been developed over a number of years with feedback and feature requests from public health laboratory epidemiolgists during its routine use in HAV outbreak surveillance.  The software aims to provide actionable results to epidemiologists with the least amount of analyst interaction despite inherent imperfections in HAV sequencing data.  As the HAV genome comprises a single segment, being a positive-sense single-stranded ribonucleic acid (RNA) of only 7.5 kilobases (kb) [@RN375], analyses can be completed relatively easily, even on a low-powered desktop computer.  

The software pipeline is written in python3, implementing a number of R packages, with pipeline control managed via Rufus.  Run configurations are defined using a yaml-formatted text file.  Installation is performed using `conda`.  

# Testing and validation

Example data are packaged with `havic`.  Tests are performed by 'doing' `havic test <test_suite>`, where `<test_suite>` is any of `hav_amplicon`, `hav_wgs`, `measles_wgs` or 'hiv_amplicon'.

# Figures


# Acknowledgements


# References

