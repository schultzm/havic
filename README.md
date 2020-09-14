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

### Advance usage

### Tips and tricks

## Release history

## Frequently Asked Questions

_Can havic be used with a custom SUBJECT sequence?_

Yes.  

_Can the the subject sequence consist of multiple contigs?_

No.

<!-- _How does it scale?_

At the outset, the pipeline compiles all the input query files into a single file using `cat` (`O(n)`, i.e., linear time complexity), discards duplicate sequence IDs (`O(n)`), throws away the bad characters in the sequence headers (`O(n)`) and writes the set of sequences to file (`O(n)`).   -->
