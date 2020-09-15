# HAVIC

[![CircleCI](https://circleci.com/gh/schultzm/havic.svg?style=svg&circle-token=9d17418bb752aa29e07f95b09af106aef7cc6b02)](https://app.circleci.com/pipelines/github/schultzm/havic)

Detect **H**epatitis **A** **V**irus **I**nfection **C**lusters from virus consensus sequences.  

## Overview

## Installation

Installation requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [git](https://git-scm.com/downloads).  After installing these packages, simply do:

    git clone https://github.com/schultzm/havic.git
    cd havic
    . install.sh

The process will take up to 30 minutes and verbose output will be printed to screen during the install.

## Usage

After installing, the most basic usage of havic is to type `havic` on the command line and hit enter/return.  If the install has worked correctly, the user should see:

```{bash}
usage: havic [-h]  ...

optional arguments:
  -h, --help  show this help message and exit

Sub-commands help:
  
    detect    Detect infection clusters from cDNA or DNA consensus sequences.
    version   Print version.
    test      Run HAVIC test using pre-packaged example data.
```

### Quickstart

    havic detect
    havic version
    havic test

### General usage

### Example usage

### Advanced usage

### Tips and tricks

## Release history

Active development.  Pre-release.  
## Frequently Asked Questions

_Why the name `havic`?_

`havic` is an acronym for **H**epatitis **A** **V**irus **I**nfection **C**luster (HAVIC), the VIC doubles as an abbreviation for Victoria, Australia, where our laboratories are based.

_Who is `havic` for?_

`havic` is for molecular epidemiologists working in public health laboratories who want to discover infection clusters in their virus sample cDNA or DNA sequences.  

_What is `havic` for?_

`havic` is for bioinformatic analysis of Hepatitis A Virus genome sequences.  It takes fasta files as input (QUERIES), maps the QUERIES to a reference (SUBJECT), extracts the alignment from the binary alignment map (bam) file, infers a phylogenetic tree from the alignment, picks infection clusters within the QUERIES using the tree and alignment as evidence.  Theoretically, `havic` can be used on other viral genomes though testing on non-HAV samples has so far been limited to Measles and SARS-CoV-2.

_How do you define SUBJECT and QUERY sequences?

To maintain consistency with already established methods, SUBJECT ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=References) nomenclature) is used interchangeably with REFERENCE, REF or reference allele [.vcf standard](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/).  SUBJECT is the backbone onto which all QUERY sequences will be mapped.  QUERY (BLAST nomenclature) is used interchangeably with ALTERNATE or ALT or alternate allele (.vcf standard).  

_Can havic be used with a custom SUBJECT sequence?_

Yes.  The HAVIC pipeline is expected to work for any non-segmented virus genome.

_Can the the SUBJECT file consist of multiple contigs?_

No.  The SUBJECT sequence needs to be a single consensus sequence from a single sample.  

_Can input QUERY samples be comprised of multiple consensus sequences from the same sample?_

No.  A QUERY file may NOT consist of multiple contigs from the same sample.  However, a QUERY file may consist of multiple sequences, one sequence from each sample.  

_Can input QUERY files consist of multiple sequences?_

Yes.  A QUERY file may either be a single consensus sequence from a single sample, or multiple samples with a single consensus sequence for each sample.  A single QUERY file can be input to `havic`, but the program is designed to accept as many QUERY files as you wish to feed it.  

_What's all this talk about consensus sequences?  I'm used to talking about contigs.

In the 2020 pandemic era, virus genome sequencing is dominated by tiled-PCR-amplicon Illumina paired-end sequencing and/or Oxford Nanpore Technologies (ONT) long read sequencing.  The typically low input nucleic acid quantity from clinical samples means that Illumina sequencing of tiled PCR amplicons is the preferred method whole genome sequencing of clinical virus samples.  Tiled amplicon Illumina sequencing allows mapping of reads from a single sample to a single reference, with the final sample genome sequence called as the consensus variants against the reference, padded by inter-variant reference bases.  The final sample sequence is not produced from a de novo assembly of reads so is referred to as a consensus sequence.  Further, in diagnostic laboratories worldwide, quantitative Reverse Transcriptase Real-time PCR (qRT-PCR, qPCR or sometimes just RT-PCR) is used to detect positive cases.  Due to difficulties associated with whole genome sequencing, diagnostic laboratorie often use Sanger sequencing of PCR products to call the strain of virus.  `havic` was originally written to discover and characterise outbreak clusters from short amplicon Sanger sequences, but now is also capable of analysis virus whole genome consensus sequences.  

_Will HAVIC work on organisms other than viruses?

Probably.  HAVIC has been designed and tested specifically to work on Hepatitis A Virus (HAV) genomes.  However, `havic` should work on any non-segmented virus genome, and successful test analyses have been performed on Measles and SARS-CoV-2 genomes.  Ultimately it is up to the analyst to decide whether `havic`'s treatment of the data makes biological sense.  
