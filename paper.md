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
  - name: Karolina Mercoula
    affiliation: 1
  - name: Michelle Sait
    orcid: 0000-0001-8782-6128
    affiliation: 1
  - name: Susan A. Ballard
    orcid: 0000-0002-0096-9474
    affiliation: 1
  - Lilly Yuen
    orcid: 0000-0003-2934-0167
    affiliation: 2
  - Leon Caly
    orcid: 0000-0003-2846-6729
    affiliation: 2
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

Hepatitis A is a viral liver disease associated with high morbidity.  Though the mortality rate is low, infections are costly and painful oftentimes resulting in hospitalisation during the acute phase of infection.  Tens of millions of people are infected with Hepatitis A Virus (HAV) each year [@RN406].  

HAV is transmitted via the faecal-oral pathway, either directly from person-to-person or through contaminated food and water [@RN398].  The virus has its highest prevalence in low-income countries; however, sporadic outbreaks occur also in high-income countries [@RN406].  Molecular epidemiological analyses of outbreak genome sequence data are essential for public health interventions aimed at attenuating outbreaks [@RN407, @RN410].  

# Statement of need

Hepatovirus A belongs to the family Picornaviridae.  The genome is a positive-sense single-stranded ribonucleic acid (RNA) approximately 7.5 kilobases (kb) in length [@RN375].  Genetic diversity of HAV is well characterised, with one serotype and seven genotypes (I-VII); of the latter, four (I, II, III and VII) infect humans, and these can be further divided into six subtypes (IA, IB, II, IIIA, IIIB, VII) [@RN404, @RN403].  

The gold standard HAV genotyping protocol [see HAVNET @RN316] has been to sequence a short 460 nucleotide (nt) amplicon spanning the junction of the VP1 and P2A regions.  With this method, sequencing is typically performed using the Sanger method with a consensus target sequence inferred from a stack of one or more overlapping sub-sequences.  In recent years, however, laboratory methods have applied variations to the HAVNET protocol [e.g., see @RN408] and whole genome sequencing will likely become the new standard [e.g., see @RN407].  

Variations in laboratory methods mean a consensus sequence may not always co-locate with the preferred genome target, especially when sequences are collected from databases compiled over years from many laboratories [e.g., see @RN411] [add plot here].  Additionally, artefactual nucleotide variation (e.g., low quality, indels, unwitting alignment of reverse and forward complement sequences) may be present.  Pairwise sequence comparisons to infer relatedness of samples is inadequate to acknowledge these uncertainties but pairwise genetic distance based inferences are often the mode [e.g., @RN407].  Further, Comparisons Hence, there is a need for tools capable of handling this variation which additionally 

Hence, there is a  or  .   tend to be made For HAV infection cluster detection, this variation presents as  


# Statement of need


The alignment of partial HAV genome sequences with the ultimate aim of detecting outbreak clusters is often troublesome.  Frame-shift mutations, incorrect orientation to the reference, and difficulty locating the true open reading frame means molecular epidemiology of HAV has hitherto been performed using the closed-source HAVNET genotyping service, NCBI BLAST [@RN172] on a case-by-case basis.  There is a need for tools that allow open-source, high-throughput, objective identification of outbreak infection clusters.

# havic

`havic` was written out of the need 

# Validation data

[@RN407, @RN408]