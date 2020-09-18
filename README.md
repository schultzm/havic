# havic

[![CircleCI](https://circleci.com/gh/schultzm/havic.svg?style=svg&circle-token=9d17418bb752aa29e07f95b09af106aef7cc6b02)](https://app.circleci.com/pipelines/github/schultzm/havic)

Detect **H**epatitis **A** **V**irus **I**nfection **C**lusters from virus consensus sequences.  

## Overview

`havic` is a bioinformatics pipeline for detecting infection clusters in Hepatitis A Virus samples from DNA or cDNA sequence data.  The pipeline is written in `python3` and uses `ruffus` to connect a number of open-source software tools to achieve this task.  The user feeds `havic` some query files via a `yaml` config file, waits for the program to run and then checks the output folder for results.  The figure below is a schematic representation of the pipeline.  

![Pipeline](https://github.com/schultzm/havic/blob/docs/havic/data/pipeline_graph.svg?raw=true)

The above pipeline is summarised briefly here.  Firstly an output directory is created to receive the output files from a `havic` run.  Query sequences are collected into a single set, duplicates sequences are discarded (based on the sequence headers in the query fasta files) and 'troublesome' characters in sequence headers are replaced with underscore.  `havic` was originally designed for analysis of the the VP1/P2A amplicon, which is the genomic marker recommended by the Hepatitis A Virus Network ([HAVNET](https://www.rivm.nl/en/havnet)).  The VP1/P2A amplicon target is the product of a nested PCR reaction, and is shown here in the context of the HAV genome:

![Amplicon](https://github.com/schultzm/havic/blob/docs/havic/data/VP1P2A.png?raw=true "The HAV genome with HAVNET amplicon, sourced from RIVM")

## Installation

Installation of `havic` requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [git](https://git-scm.com/downloads).  After installing these packages, simply do:

    git clone https://github.com/schultzm/havic.git
    cd havic
    . install.sh

The installation process will take up to 30 minutes with verbose output printed to screen during the install.  If the installation fails, read the screen output to determine the error via traceback.  Submit installation issues to github.  Installation has been tested via continuous integration on CircleCI and tested inside a conda environment.  

## Usage

### Quickstart

After installing, activate the conda environment by doing `conda activate havic_env`.  The most basic usage of `havic` is to type `havic` on the command line and hit enter/return.  If the install has worked correctly, the user should see:

    usage: havic [-h]  ...

    optional arguments:
    -h, --help  show this help message and exit

    Sub-commands help:
    
        detect    Detect infection clusters from cDNA or DNA consensus sequences.
        version   Print version.
        test      Run havic test using pre-packaged example data.

The program is accessed via three subcommands, with help via the `-h` suffix.  

`havic detect` is the main sub-command.  Use this for detecting infection clusters from user-specified cDNA or DNA consensus sequences.  
`havic version` will print the installed version to `stdout`.  
`havic test` will run `havic detect` on a pre-packaged test dataset.  If successful, the analyst should see `ok` at the end of each test.

### Example usage

For a basic analysis of HAV VP1/P2A amplicons, the analyst will likely need to view the output in the context of circulating HAV strains.  Hence, the first step of analysis should be to collect samples from a database of sequences (e.g., NCBI GenBank or HAVNET).  `havic` requires at least three query sequences to run.  After collecting the sequences into a fasta file/s, edit the `yaml` config file before trying to perform any analysis.

#### Editing the `yaml` file for parsing by `havic detect`

`havic detect` receives instructions from a `yaml` config file via the command `havic detect path/to/yaml.yaml`.  The `test.yaml` file from `havic/havic/data/havic_detect.yaml` is presented below as an example:

    ---
    FORCE_OVERWRITE_AND_RE_RUN:
    Yes # Yes for full re-run, No to start from an interrupted run,

    DEFAULT_REFS:
    Yes # Yes if using havic pre-packaged SUBJECT test data, No otherwise

    DEFAULT_QUERIES:
    Yes # Yes if using havic pre-packaged QUERY test data, No otherwise

    SUBJECT_FILE: # the "SUBJECT" sequence in BLAST terms, i.e., reference genome
    data/NC_001489.fa # relative or absolute paths to fasta file
    # if DEFAULT_REFS is Yes, path will be prefixed to use pre-packaged data

    SUBJECT_TARGET_REGION: # the target region of the genome to focus on
    data/havnet_amplicon.fa # in fasta format, relative or absolute paths okay
    # if DEFAULT_REFS is Yes, path will be prefixed to use pre-packaged data

    OUTDIR: # the parent directory for the results folders
    havic_test_results/r1 # relative or absolute path to parent result folder

    TREE_ROOT:
    midpoint # sequence name to root iqtree on, or midpoint for midpoint root

    RUN_PREFIX:
    HAV_all_

    PLOTS:
    Yes # Yes to make plots (slow for large runs), No otherwise.

    CLUSTER_PICKER_SETTINGS: # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4228337/
    executable:
        ClusterPicker
    coarse_subtree_support: # divide tree into subtrees at/above this threshold
        70
    fine_cluster_support: # branch support minimum value for clusters of tips
        95
    distance_fraction: # float please, genetic distance
        0.01 # (e.g., 1 SNP in 100 bp = 0.01)
    large_cluster_threshold:
        15
    distance_method:
        valid # options are ambiguity, valid, gap, or abs

    IQTREE2_SETTINGS:
    executable:
        iqtree # command to call iqtree2
    threads: # threads
        '-T AUTO -ntmax 24' # automatically determine the best threading parameter
    model_finder: # model-finder
        '-m MFP' # extended model find with FreeRate heterogeneity + tree inference
    state_frequency:
        '+FO' # Optimized sgate frequencies by maximum-likelihood
    ultrafast_bootstrap:
        '--ufboot 1000' # use of aLRT will cause ClusterPicker to fall over
    protect_violations: #to protect against severe model violations
        '--bnni'
    redo: # recompute everything in iqtree run
        '--redo' # leave empty string ('') if not wanting to redo, else '--redo'

    MAPPER_SETTINGS:
    executable:
        minimap2
    other:
        --secondary=no -Y
    k_mer: # select an odd number, between 3 and 27 inclusive
        -k 5 # 5 has been good for the HAV amplicon seqs, adjust sensibly

    HIGHLIGHT_TIP:
    - CmvAXJTIqH # Specify tip name to highlight in final plot
    - CCHkiFhcxG # Specify tip name to highlight in final plot
    - PAvYXhYkLM # Specify tip name to highlight in final plot

    TRIM_SEQS: # these sequences will be trimmed to length of SUBJECT_AMPLICON
    - AY644337_55443_seq_1 # these are sequences in the QUERY_FILES
    - RIVM-HAV171_64913_seq_2_MapsOutsideTrimRegionSoEmpty
    - nDNLdjtgha#HashInSeqName
    - '' # give it nothing
    - xyzyx # give it a non-name

    QUERY_FILES:
    - data/example1.fa # relative or absolute paths to fasta files
    - data/example2.fa
    - xyz # to test a dud file name
    - '' # to test an empty file name (which would return a folder, not file)
    ...

Before starting a run, `cd` to a working directory (preferably not inside the git cloned folder).  Either copy the above `yaml` to file, or use `wget https://raw.githubusercontent.com/schultzm/havic/master/havic/data/havic_detect.yaml`.  For more information on the `yaml` standard, refer to [https://yaml.org/](https://yaml.org/).  

Lets go through the `yaml` step-by-step.

##### Opening and closing fields, nesting

`yaml` code blocks open and close with `---` and `...`, respectively.  Ensure your file includes these lines.  Indents are two spaces to increase a nesting level.  

##### Force overwrite and re-run

    FORCE_OVERWRITE_AND_RE_RUN:
    Yes

`havic` manages tasks via [`ruffus`](https://code.google.com/archive/p/ruffus/), and out-of-date stages of the pipeline will be re-run as required.  To start a new run or force overwrite files in the OUTDIR, set `FORCE_OVERWRITE_AND_RE_RUN` to `Yes`.  Otherwise to start off from the last point, set to `No`.  

##### Default Subject and Queries

    DEFAULT_SUBJECT:
    Yes # Yes if using havic pre-packaged SUBJECT (i.e., 'reference') sequence and region test data, No otherwise

    DEFAULT_QUERIES:
    Yes # Yes if using havic pre-packaged QUERY test data, No otherwise

If `Yes`, for `DEFAULT_SUBJECT`, `havic` will prefix the filepaths in `SUBJECT_FILE` `SUBJECT_TARGET_REGION` with the `havic` install path (using `pkg_resources.resource_filename`) object to search in the installation directory for the pre-packaged data.  If `DEFAULT_SUBJECT` is set to `No`, then the filepaths provided .  If `DEFAULT_QUERY` is set to `No`, then 

##### Subject/Reference sequence

    SUBJECT_FILE: # the "SUBJECT" sequence in BLAST terms, i.e., reference genome
    data/NC_001489.fa # relative or absolute paths to fasta file
    # if DEFAULT_SUBJECT is Yes, path will be prefixed to use pre-packaged data

`havic` will use this fasta sequence as the subject/reference sequence.  If a different reference is required, change the path value.  

##### Subject target region

    SUBJECT_TARGET_REGION: # the target region of the genome to focus on
    data/havnet_amplicon.fa # in fasta format, relative or absolute paths okay
    # if DEFAULT_REFS is Yes, path will be prefixed to use pre-packaged data

This regions will guide trimming of the alignment.  In this example, the VP1/P2A region is the target region.  Sample names listed in TRIM_SEQS will be trimmed to match the boundaries of this region.  A sequence is used here instead of a bed coordinates file because the exact boundaries of the target region in the final alignment are not always obvious.  After mapping this region to the subject sequence, the boundaries become obvious.  Automatic delineation of this region alleviates the need for the analyst to manually search for and define the boundaries.  

##### Output directory

    OUTDIR: # the parent directory for the results folders
    havic_test_results/r1 # relative or absolute path to parent result folder

Specify the path to the output directory.  The files listed in the table below will be sent to this directory as the run progresses.  

###### Output files

Stage number | Stage name | File or directory name
---:|:---|:---
1 | create_outdir | `havic_test_results/r1`
2 | compile_input_fasta | `HAV_all_duplicate_seqs.txt`
2 | compile_input_fasta | `HAV_all_seq_id_replace.tsv`
2 | compile_input_fasta | `HAV_all_tmpfasta.fa`
3 | map_input_fasta_to_ref | `HAV_all_map.bam`
3 | map_input_fasta_to_ref | `HAV_all_map.bam.bai`
4 | bam2fasta | `HAV_all_map.bam2fasta.R`
4 | bam2fasta | `HAV_all_map.bam2fasta.Rout`
4 | bam2fasta | `HAV_all_map.stack.fa`
5 | get_cleaned_fasta | `HAV_all_map.stack.trimmed.fa`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.bionj`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.ckp.gz`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.contree`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.iqtree`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.log`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.mldist`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.model.gz`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.splits.nex`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.treefile`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.ufboot`
6 | run_iqtree | `HAV_all_map.stack.trimmed.fa.uniqueseq.phy`
7 | root_iqtree | `HAV_all_map.stack.trimmed.fa.rooted.treefile`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa_HAV_all_map.stack.trimmed.fa.rooted_clusterPicks_cluster4_sequenceList.txt`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa_HAV_all_map.stack.trimmed.fa.rooted_clusterPicks.fas`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa.rooted_clusterPicks_list.txt`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa.rooted_clusterPicks_log.txt`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa.rooted_clusterPicks.nwk`
8 | clusterpick_from_rooted_iqtree_and_cleaned_fasta | `HAV_all_map.stack.trimmed.fa.rooted_clusterPicks.nwk.figTree`
9 | summarise_cluster_assignments | `HAV_all_map.stack.trimmed.fa.rooted_clusterPicks_summarised.txt`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa_SNPcountsOverAlignLength.csv`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa_SNPdists.csv`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa_SNPdists.pdf`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa.rooted.treefile_1percent_divergence_valid_msa.pdf`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa.Rplot.R`
10 | plot_results_ggtree | `HAV_all_map.stack.trimmed.fa.Rplot.Rout`
11 | pipeline_printout_graph | `pipeline_graph.svg`

##### Setting the location of the tree root

    TREE_ROOT:
    midpoint # sequence name to root iqtree on, or midpoint for midpoint root

For visual representation only, the tree root is set to orientate the plot in `HAV_all_map.stack.trimmed.fa.rooted.treefile_1percent_divergence_valid_msa.pdf`.  The tree root does not affect cluster definitions.  

##### Set the prefix of output filenames

    RUN_PREFIX:
    HAV_all_

To facilitate tracking of output files, the user is able to specify a custom prefix for output files.  

##### Draw results plots

    PLOTS:
    Yes # Yes to make plots (slow for large runs), No otherwise.

This setting controls the drawing of output plots.  The plots are helpful to understand how the multiple sequence alignment affects tree topology, cluster detection and pairwise SNP distances.

![Heatmap](https://github.com/schultzm/havic/blob/docs/havic/data/heatmap_SNPs?raw=true)

![Tree](https://github.com/schultzm/havic/blob/docs/havic/data/tree_MSA_clusters.png?raw=true)


##### Input query files

Input query sequences should be in fasta format with one sequence per sample.  Multiple samples may be included per file, and/or multiple files may be passed to `havic`.  Query sequences within files will be reverse complemented as necessary during their mapping to the subject/reference.  If the query sequence files are named `batch1.fa`, `batch2.fa`, `batch3.fa`,  edit the `QUERY_FILES` section of the `yaml` file as follows:

    QUERY_FILES:
    - batch1.fa # relative or absolute paths to fasta files
    - batch2.fa
    - batch3.fa

##### Highlighting samples of interest

To highlight query sequences in the final plots, list the sequence names under `HIGHLIGHT_TIP` in the `yaml`, otherwise ignore this section.  

##### Trimming sequences to genomic region of interest

To trim input queries to the reference VP1/P2A amplicon, list the sequence name of the query under `TRIM_SEQS`, otherwise ignore this section.  


### Advanced usage

preifx with something != RUN_PREFIX
 During development of `havic`, it was recognised that HAV surveillance will likely move to whole genome sequencing in the near future.  To improve utility of `havic` over the coming years, the software is written to allow the user to pass in any query sequence and any subject sequence.  Prior to phylogenetic analysis, query headers listed under `TRIM_SEQS` will be trimmed to the subject target region given by `SUBJECT_TARGET_REGION`.  Defining the subject target region in this way as opposed to using a bed file of coordinates is a feature to allow the user to not have to know in advance exactly where the target region is.  With changing references, the target region may differ slightly from the expected region, so by allowing the mapper to find the region, the user is relieved the burden of having to find and define the region a priori.  

`havic` will read query sequences and map them to the subject sequence provided under `SUBJECT_FILE`.  The subject file can be any single contig the user desires.  `havic` has been tested on HAV (~7.5kb) genomes with preliminary testing also being successful for Measles (~15.9kb) and SARS-CoV-2 (~30kb) genomes.  The upper limit has not yet been found.  


### Tips and tricks

## Release history

Active development.  Pre-release.  
## Frequently Asked Questions

_Why the name_ `havic`_?_

`havic` is an acronym for **H**epatitis **A** **V**irus **I**nfection **C**luster (HAVIC), the **VIC** acknowledges that the development team hails from Victoria, Australia.

_Who is `havic` for?_

`havic` is for molecular epidemiologists working in public health laboratories who want to discover infection clusters in their virus sample cDNA or DNA sequences.  

_What is `havic` for?_

`havic` is for bioinformatic analysis of Hepatitis A Virus genome sequences.  It takes fasta files as input (QUERIES), maps the QUERIES to a reference (SUBJECT), extracts the alignment from the binary alignment map (bam) file, infers a phylogenetic tree from the alignment, picks infection clusters within the QUERIES using the tree and alignment as evidence.  Theoretically, `havic` can be used on other viral genomes though testing on non-HAV samples has so far been limited to Measles and SARS-CoV-2.

_How do you define SUBJECT and QUERY sequences?_

To maintain consistency with already established methods, SUBJECT ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=References) nomenclature) is used interchangeably with REFERENCE, REF or reference allele [.vcf standard](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/).  SUBJECT is the backbone onto which all QUERY sequences will be mapped.  QUERY (BLAST nomenclature) is used interchangeably with ALTERNATE or ALT or alternate allele (.vcf standard).  

_Can havic be used with a custom SUBJECT sequence?_

Yes.  The havic pipeline is expected to work for any non-segmented virus genome.

_Can the the SUBJECT file consist of multiple contigs?_

No.  The SUBJECT sequence needs to be a single consensus sequence from a single sample.  

_Can input QUERY samples be comprised of multiple consensus sequences from the same sample?_

No.  A QUERY file may NOT consist of multiple contigs from the same sample.  However, a QUERY file may consist of multiple sequences, one sequence from each sample.  

_Can input QUERY files consist of multiple sequences?_

Yes.  A QUERY file may either be a single consensus sequence from a single sample, or multiple samples with a single consensus sequence for each sample.  A single QUERY file can be input to `havic`, but the program is designed to accept as many QUERY files as you wish to feed it.  

_What's all this talk about consensus sequences?  I'm used to talking about contigs._

In the 2020 pandemic era, virus genome sequencing is dominated by tiled-PCR-amplicon Illumina paired-end sequencing and/or Oxford Nanpore Technologies (ONT) long read sequencing.  The typically low input nucleic acid quantity from clinical samples means that Illumina sequencing of tiled PCR amplicons is the preferred method whole genome sequencing of clinical virus samples.  Tiled amplicon Illumina sequencing allows mapping of reads from a single sample to a single reference, with the final sample genome sequence called as the consensus variants against the reference, padded by inter-variant reference bases.  The final sample sequence is not produced from a de novo assembly of reads so is referred to as a consensus sequence.  Further, in diagnostic laboratories worldwide, quantitative Reverse Transcriptase Real-time PCR (qRT-PCR, qPCR or sometimes just RT-PCR) is used to detect positive cases.  Due to difficulties associated with whole genome sequencing, diagnostic laboratorie often use Sanger sequencing of PCR products to call the strain of virus.  `havic` was originally written to discover and characterise outbreak clusters from short amplicon Sanger sequences, but now is also capable of analysis virus whole genome consensus sequences.  

_Will havic work on organisms other than viruses?_

Probably.  havic has been designed and tested specifically to work on Hepatitis A Virus (HAV) genomes.  However, `havic` should work on any non-segmented virus genome, and successful test analyses have been performed on Measles and SARS-CoV-2 genomes.  Ultimately it is up to the analyst to decide whether `havic`'s treatment of the data makes biological sense.  
