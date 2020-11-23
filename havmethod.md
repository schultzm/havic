---
title: "HAVIC SOP"
author: "Mark Schultz"
date: "17/07/2020"
output: pdf_document
bibliography: references.bib
csl: harvard-the-university-of-melbourne.csl
---


#### SOP NUMBER: _**MMS168**_
#### SOP NAME: _**Hepatitis A virus infection cluster (HAVIC) detection using phylogenetic analysis of cDNA sequences**_
#### CUSTODIAN: Dr Mark B Schultz
#### AUTHORISED BY: Dr Susan Ballard
#### VERSION: 1.0

# Background & Introduction

Globally, Hepatitis A is a significant cause of morbidity.  This infectious disease of the liver is caused by the Hepatitis A Virus (HAV), which is transmitted from person-to-person and through contaminated food and water principally via the faecal-oral pathway [@RN398].  Hepatitis A is vaccine preventable and immunity to the virus develops post-infection.  In high income countries with generally higher sanitation standards, transmission of HAV among children is low so many adults may remain susceptible to infection [@RN405].  

Despite Australia being ranked as a high income country, the proportion of the population susceptible to HAV has decreased in recent decades due to high rates of vaccination and a continual influx of immigrants from HAV endemic countries [@RN395, @RN385, @RN387, @RN400].  However, the risk of occasional outbreaks remains: Australia has experienced sporadic outbreaks in the last 10 years as a result of import of contaminated food [@RN399], or via sexual transmission events [@RN401].  

Cases of HAV are nationally notifiable within Australia.  Laboratory definitive evidence is sufficient to confirm a case, which is provided by HAV nucleic acid testing [@RN402].  In attempts to standardise HAV laboratory protocols, Australian testing laboratories have sometimes adopted the gold standard amplicon sequencing protocol provided by HAVNet (the Hepatitis A Virus Network).  

HAV is a non-enveloped Hepatovirus in the family Picornaviridae.  The genome is a positive-sense single-stranded ribonucleic acid (RNA) approximately 7.5 kilobases (kb) in length [@RN375].  Genetic diversity of HAV is well characterised, with one serotype and seven genotypes (I-VII); of the latter, four (I, II, III and VII) infect humans, which can be divided further into six subtypes (IA, IB, II, IIIA, IIIB, VII) [@RN404, @RN403].  

The gold standard HAVNet genotyping protocol [@RN316] recommends sequencing of a short 460 nucleotide (nt) amplicon spanning the junction of the VP1 and P2A regions.  In short, viral nucleic acid is extracted using an RNA isolation technique.  Reverse-transcriptase (RT) enzyme is used to back transcribe the RNA to complementary deoxyribonucleic acid (cDNA).  A two-stage PCR is then applied, whereby the first stage amplifies a 503 nt fragment and a second nested PCR stage amplifies a 460 nt.  The latter fragment is DNA sequenced, typically using Sanger di-deoxy capillary sequencing.  


# Scope

The scope of this Standard Operating Procedure (SOP) is to describe the bioinformatic processes used at MDU-PHL to infer genotype and detect outbreak clusters of Hepatitis A infections from HAV VP1/P2A cDNA sequence data.  This SOP describes use of the HAV Infection Cluster (HAVIC) detection software as applied to cDNA sequence data.  Further, it describes the output of the analysis and discusses limitations and areas for improvement in the workflow.  HAVIC was written by the bioinformatics division at MDU-PHL and is freely available under open source GNU AGPL3 licence from https://github.com/schultzm/HAVIC.


# Definitions

Below are definitions for acronyms used in this document.  

| Acronym | Definition |
|---------|------------|
| AGPL3 | Affero General Public Licence version 3 |
| CLI | Command Line Interface |
| DNA | deoxyribonucleic acid |
| FOFN | File of File Names |
| HAV | Hepatitis A Virus |
| HAVIC | software for detecting Hepatitis A Virus Infection Clusters |
| IUPAC code | Refer to [web for further details](https://genome.ucsc.edu/goldenPath/help/iupac.html) |
| IUPAC | International Union of Pure and Applied Chemistry |
| MDU-PHL | Microbiological Diagnostic Unit - Public Health Laboratory |
| N | Any nucleotide (A, C, G or T) |
| NCBI | [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/) |
| PCR | Polymerase Chain Reaction |
| QHFSS | Queensland Health Forensic Science Services |
| RNA | ribonucleic acid |
| RT | Reverse Transcriptase enzyme |
| VIDRL | Victorian Infectious Diseases Reference Laboratory |
| cDNA | complementary DNA |
| kb | kilobases |
| nt | nucleotide |

# Internal documents

_**Any other documentation that this document refers to - particularly lab method or other SOP or validation documents**_

| SOP/Doc. no. | Title |
| :--- | :--- |
|**MMSXXX/VXXX** | SOP and validation documentation that has been logded with QA |

# Quality Control / Reference Material

HAVIC test suite

# Risk and Hazard Assessment

## Safety

Refer to Risk register at : E:\\groups\\ehs\\FM1187 RiskRegister Records on MDU internal server

# Safety Alert / Hazard Controls

| Hazard | Risk | Control | Associated RA\(s\) |
| :--- | :--- | :--- | :--- |
| Poor quality data analysed | Incorrect results released | Sequence data must pass QC metrics before analysis | RA197 |
| Computer workstation not ergonomically set up | Minor to serious injury |  | RA197 |
| Electrical shock | Minor to serious injury | Test and tag hardware | RA197 |

# Specimen Handling Requirements

Refer to _**MMSXXX**_ for specimen handling requirements

Samples from cases presenting with HAV symptoms are sent to either the Victorian Infectious Diseases Reference Laboratory (VIDRL) or Queensland Health Forensic and Scientific Services (QHFSS) for nucleic acid testing.  Samples are not physicaly handled at MDU-PHL.  Computer files for analysis at MDU-PHL are sent from QHFSS and VIDRL.  

# Methods

To install the HAVIC software, follow the instructions at https://github.com/schultzm/HAVIC.

## Detecting HAV infection clusters using HAVIC

A typical job request will come to the analyst via MeisterTask in a form:  

```{bash}
...could you please run the HAV pipeline. Thanks!
 
New sequences:
•	6 WA
•	2 international
```

To run the job, start a tmux session on one of the MDU-PHL servers:
```{bash}
tmux new -s hav
```

Inside the `tmux` session, run the following command:

```
seq 1 3 | while read i;
do
  echo /home/schultzm/.local/bin/havic detect -q /mnt/seq/MDU/SURVEILLANCE/HAV/{*.fa*,*.FA*} -t $(grep '>' /mnt/seq/MDU/SURVEILLANCE/HAV/*nternat* | cut -d '>' -f 2 | cut -d ' ' -f 1 | perl -pe 's/\R/ /g' | xargs) -o $(date '+%Y%m%d')/r${i} -p $(date '+%Y%m%d')_r${i}
done | parallel --bar --verbose {}
```

# Test Interpretation

_**What is the output of the test/tool and what does it mean? Include an example of the output with a bit of an explanation.**_


## Acceptance/Rejection Criteria

QC of capillary sequence chromatogram files is performed at the submitting laboratory.  A number of QC steps are performed as part of the HAVIC run.  

The reference amplicon spans positions 2915 to 3374 on NCBI genome `NC_001489.1`.  The reference genotype is IB.  The `.bed` formatted amplicon is:  

```{bash}
NC_001489_1	2914	3374
```

The reference amplicon is stored in `HAVIC/havic/data/havnet_amplicon.py` and the reference genome is stored in `HAVIC/havic/data/NC_001489.fa`.  If the reference bed file is `ref.bed`, then the command to extract the reference amplicon from the full reference sequence is:  

```{bash}
bedtools getfasta -fi NC_001489.fa -bed ref.bed
```

Sequences with FASTA-format non-desirable characters are filtered as follows:  
- replace any non `A-Z`, `a-z`, `0-9` with `_`  
- remove `_(reversed)`  

Duplicate sequences are discarded with duplicate sequences infered from fasta header (as opposed to cDNA sequence identity).  

## Internal Quality Control

To test the HAVIC software installation do `havic test -c $(which clusterpicker)`.  This will run the software on the test dataset.  

To troubleshoot and check HAVIC software dependencies do `havic depcheck`.  

# Troubleshooting

| Problem | Solution |
| :--- | :--- |
| The list of input files is too long for passing in via CLI | Use a File of File Names (FOFN) (this is not yet implemented) |
| The list of fasta IDs to trim is too long for passing in via CLI | Use a file of trim IDs |  
| ClusterPicker results are non-sensical | - re-run with adjusted ClusterPicker settings |
| ClusterPicker results are non-sensical | - swap out ClusterPicker for Phydelity |
| ClusterPicker is not running due to `Error: Invalid or corrupt jarfile ClusterPicker.jar` | Point havic to the correct ClusterPicker.jar file using `havic detect -c path/to/ClusterPicker.jar` |
| Miscellaneous | Submit issue to https://github.com/schultzm/HAVIC/issues |
| HAVIC crashes prematurely | Check HAVIC log files for traceback |
| Unable to find input sequence in output files | Check log files.  Output sequence has likely had undesirable characters replaced with `_` |  
| HAVIC takes a very long time to run | Reduce the number of input sequences.  For 2000 input sequences, runs can take 8 hours or more |
| There are no graphical output files in the output folder but I used the `-m` option | Check the `*.Rout` files for traceback |


# Reporting

At the completion of three independent runs, results are reported back to the requestor via the original MeisterTask request.  This is done by commenting in the request with the location of the results folders (e.g., "results are now available at `destination/dir/r{1..3}`").  After commenting, move the task into `Done` and assign the task to the nominated recipient.

# Record Management

* Refer to BIS001 (Sequence QC) for details on record management of bioinformatics data files.  
* input files are stored in `/home/seq/MDU/SURVEILLANCE/HAV/{*.fa*,*.FA*}`
* international match files are stored in `/home/seq/MDU/SURVEILLANCE/HAV/*nternational*`

# Quality Parameters

## Training

* Assessment of knowledge by reading this SOP  
* UNIX CLI competency required  
* Completion of `Introduction to MDU Server` induction  
* Completion of `Running MDU jobs` induction  
* Evidence of successfully running this tool under supervision  


## Test Performance Characteristics

Refer to VXX Validation and Verification Protocol for SOP MMSXXX.

## Measurement Uncertainty

_**This is normally only required for quantitative tests**_

## Proficiency Testing Program

_**Check these with MDDI**_

# Document History

| Version number	 | Prepared/Review by | Date | Brief description of changes made to this version | Re-verification/validation performed Yes/No/NA | Authorised by | 
|-------|-------|-----|-----|------|-----|
| Draft | Dr Mark B Schultz | 2020-07-22 | First version | NA | Director |


# References
 
