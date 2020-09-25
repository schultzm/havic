---
title: 'havic: detection of Hepatitis A Virus Infection Clusters from clinical cDNA sequences'
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
  - name: Patiyan Andersson
    orcid: 0000-0003-2705-1847
    affiliation: 1
  - name: Marion Easton
    affiliation: 3
  - name: Linda T. Viberg
    orcid: 0000-0002-4174-1142
    affiliation: 1
  - name: William Pitchers
    orcid: 0000-0003-0385-5939
    affiliation: 1
  - name: Karolina Mercoulia ***
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
  - Mike Catton
  - Julian Druce
  - Courtney Lane
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

Tens of millions of people are infected with Hepatitis A Virus (HAV) each year [@RN406].  Though the mortality rate is low, this viral liver disease is associated with high morbidity.  The economic burden of the disease is high, as cases often required hospitalisation during the acute phase of infection.  

HAV is transmitted via the faecal-oral pathway, either directly from person-to-person or through contaminated food and water [@RN398].  The virus has its highest prevalence in low-income countries; however, sporadic outbreaks do occur in high-income countries [@RN406] typically via import of frozen foods from HAV endemic countries [@RN410].  Molecular epidemiological analyses of outbreak genome sequences are essential for public health interventions aimed at attenuating outbreaks [@RN407, @RN410].  `havic` is a tool written for molecular epidemiological purposes and has been developed for routine public health molecular epidemiological analysis and outbreak surveillance.  The software has been tested over a number of years and iterations have added feature requests from front-line public health epidemiologists involved in contact tracing and case history investigations.  

# Statement of need

Hepatovirus A belongs to the family Picornaviridae.  The genome is a positive-sense single-stranded ribonucleic acid (RNA) approximately 7.5 kilobases (kb) in length [@RN375].  Genetic diversity of HAV is well characterised, with one serotype and seven genotypes (I-VII); of the latter, four (I, II, III and VII) infect humans, and these can be further divided into six subtypes (IA, IB, II, IIIA, IIIB, VII) [@RN404, @RN403].  

The gold standard genotyping protocol for HAV cases [see HAVNET @RN316] has been to sequence a short 460 nucleotide (nt) amplicon spanning the VP1 and P2A junction.  Sequencing is typically performed using the Sanger method forming a consensus target sequence from a stack of one or more overlapping sub-sequences.  Pairwise comparisons are then made to public databases on NCBI or HAVNET to infer relatedness of samples and strain origins [e.g., see @RN407 and @RN410].    

In recent years, however, laboratory methods have applied variations to the HAVNET protocol [e.g., see @RN408] and whole genome sequencing is becoming the new standard [e.g., see @RN407].  

Variations in laboratory methods mean a consensus sequence may not always co-locate within the same genome target, especially when sequences are collected from databases compiled over many years from myriad laboratories [e.g., see @RN411] [add plot here].  Additionally, artefactual nucleotide variation (e.g., low quality, false indels, incompatible orientation) may be present.  Pairwise sequence comparisons to infer relatedness of samples from genetic distance is often the mode for recovering epidemiological links ; however, crude distance based methods are outdated.  Hence, there is a need for bioinformatic tools that can cope with these challenges.

`havic` was written by bionformaticists from a public health laboratory who were tasked with the routine analysis of HAV sequence data for epidemiological purposes.  The software aims to provide actionable results to epidemiologists with the least amount of fuss despite the inherent imperfections in HAV sequence data.  


# Validation data

Example data is packaged with the software.  The user can run the tests by doing `havic test <test_suite>`, where `<test_suite>` can be `hav_amplicon`, `hav_wgs` or `measles_wgs`.  was collected [@RN407, @RN408]

# Figures


# Acknowledgements


# References

