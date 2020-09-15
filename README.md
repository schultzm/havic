# HAVIC

[![CircleCI](https://circleci.com/gh/schultzm/havic.svg?style=svg&circle-token=9d17418bb752aa29e07f95b09af106aef7cc6b02)](https://app.circleci.com/pipelines/github/schultzm/havic)

Detect **H**epatitis **A** **V**irus **I**nfection **C**lusters from HAVNET amplicon sequences.  

## Overview

## Installation

Installation requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [git](https://git-scm.com/downloads).  After installing these packages, simply do:

    git clone https://github.com/schultzm/havic.git
    cd havic
    . install.sh

The process will take up to 30 minutes and you should see verbose output during the install.

## Usage

### Quickstart

    havic detect
    havic version
    havic test

### General usage

### Example usage

### Advanced usage

### Tips and tricks

## Release history

## Frequently Asked Questions

_What is `havic` for?_

`havic` is for bioinformatic analysis of Hepatitis A Virus genomes. 

_How do you define SUBJECT and QUERY sequences?

To maintain consistency with already established methods, SUBJECT ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=References) nomenclature) is used interchangeably with REFERENCE, REF or reference allele [.vcf standard](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/).  SUBJECT is the backbone onto which all QUERY sequences will be mapped.  QUERY (BLAST nomenclature) is used interchangeably with ALTERNATE or ALT or alternate allele (.vcf standard).  

_Can havic be used with a custom SUBJECT sequence?_

Yes.  The HAVIC pipeline is expected to work for any non-segmented virus genome.

_Can the the SUBJECT file consist of multiple contigs?_

No.  The SUBJECT sequence needs to be a single contig from a single sample.  That is, for any analysis, the SUBJECT s

_Can input QUERY files consist of multiple contigs?_

No.  A QUERY file may NOT consist of multiple contigs from the same sample.  However, a QUERY file may consist of multiple sequences, one from each sample.  

_Can input QUERY files consist of multiple sequences?_

Yes.  A QUERY file may either be a single contig from a single sample, or multiple samples with a single contig for each sample.  

_Will HAVIC work on organisms other than viruses?

HAVIC has been designed and tested to work on Hepatitis A Virus (HAV) genomes.  Theoretically `havic` should work on any non-segmented virus genome.  `havic` has succeeded in test analyses of Measles and SARS-CoV-2 genomes.  

<!-- _How does it scale?_

At the outset, the pipeline compiles all the input query files into a single file using `cat` (`O(n)`, i.e., linear time complexity), discards duplicate sequence IDs (`O(n)`), throws away the bad characters in the sequence headers (`O(n)`) and writes the set of sequences to file (`O(n)`).   -->
