---
title: 'havic: detect Hepatitis A Virus Infection Clusters in clinical sample sequences'
tags:
  - Python3
  - Python
  - bioinformatics
  - pipeline
  - ruffus
  - mapping
  - consensus
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

Globally, Hepatitis A Virus (HAV) infections are a significant cause of morbidity.  HAV is transmitted via the faecal-oral pathway, from person-to-person and through contaminated food and water [@RN398].  HAV is prevalent in countries with poor sanitation.  nations with poor sanitation  HAV is a non-enveloped Hepatovirus in the family Picornaviridae.  The genome is a positive-sense single-stranded ribonucleic acid (RNA) approximately 7.5 kilobases (kb) in length [@RN375].  Genotyping Genetic diversity of HAV is well characterised, with one serotype and seven genotypes (I-VII); of the latter, four (I, II, III and VII) infect humans, which can be divided further into six subtypes (IA, IB, II, IIIA, IIIB, VII) [@RN404, @RN403].  

The gold standard HAVNET genotyping protocol [@RN316] recommends sequencing of a short 460 nucleotide (nt) amplicon spanning the junction of the VP1 and P2A regions.  In short, viral nucleic acid is extracted using an RNA isolation technique.  Reverse-transcriptase (RT) enzyme is used to back transcribe the RNA to complementary deoxyribonucleic acid (cDNA).  A two-stage PCR is then applied, whereby the first stage amplifies a 503 nt fragment and a second nested PCR stage amplifies a 460 nt.  The latter fragment is DNA sequenced, typically using Sanger di-deoxy capillary sequencing.  Hepatitis A is vaccine-preventable and immunity to the virus develops post-infection.  In high income countries with generally higher sanitation standards, transmission of HAV among children is low so many adults may remain susceptible to infection [@RN405].  

Despite Australia being ranked as a high income country, the proportion of the population susceptible to HAV has decreased in recent decades due to high rates of vaccination and a continual influx of immigrants from HAV endemic countries [@RN395, @RN385, @RN387, @RN400].  However, the risk of occasional outbreaks remains: Australia has experienced sporadic outbreaks in the last 10 years as a result of import of contaminated food [@RN399], or via sexual transmission events [@RN401].  



Cases of HAV are nationally notifiable within Australia.  Laboratory definitive evidence is sufficient to confirm a case, which is provided by HAV nucleic acid testing [@RN402].  In attempts to standardise HAV laboratory protocols, Australian testing laboratories have sometimes adopted the gold standard amplicon sequencing protocol provided by HAVNET (the Hepatitis A Virus Network).  
