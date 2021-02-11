# apa: alternative polyadenylation (APA) analysis

* [About](#about)
* [Run with Docker](#run-with-docker)
  * [Build the Docker image](#build-the-docker-image)
  * [Example](#example)
    * [Data folder organization and structure](#data-folder-organization-and-structure)
      * [Structure of annotation.tab](#structure-of-annotation.tab)
      * [Structure of example.config](#structure-of-example.config)
    * [Running the Example](#running-the-example)
* [Installation as standalone](#installation-as-standalone)
* [Documentation](#documentation)
  * [Download and prepare example library](#download-and-prepare-example-library)
  * [Process example library](#process-example-library)
  * [Prepare and process example control vs test comparison](#prepare-and-process-example-control-vs-test-comparison)
* [Authors](#authors)
* [Reporting problems](#reporting-problems)

## About

apa is a Python framework for processing and analysing 3'-end targeted sequence data to study alternative polyadenylation. [apa](https://github.com/grexor/apa) depends on [pybio](https://github.com/grexor/pybio) (basic handling of annotated genomes), [RNAmotifs2](https://github.com/grexor/rnamotifs2) (analysis of regulatory motif clusters) and other open-source software (DEXSeq, STAR short-read aligner).

The inclusive nature of the framework, together with novel integrative solutions (differential polyA site usage and RNA-protein binding via RNA-maps, cluster motif analysis), results in the following computational capabilities:

+ management of diverse high-throughput sequencing datasets (pre-processing, alignment, annotation),
+ polyA site database (atlas) construction and comparison to existing polyA resources,
+ identification of genes that undergo alternative polyadenylation (DEXSeq),
+ identification of motifs influencing polyA site choice (RNAmotifs2),
+ identification of motifs influencing alternative splicing (DEXSeq and RNAmotifs2),
+ integration with iCLIP (RNA-protein binding) and computing RNA-maps,
+ and other.

## Run with Docker

Since this python package has several dependencies (`pybio` for genome download and manupulation, for example), the easiest way to try out the package with an example dataset is by running the docker image (provided with the [Dockerfile](Dockerfile)).

### Build the Docker image

Clone this repository (`git clone https://github.com/grexor/apa.git`) and run `build.sh` to build the Docker image (you need to have Docker installed).

This will build a Docker image with *apa*, *pybio* and all other dependencies installed. It will also create a "data" folder on your local drive (inside the same folder where `build.sh` is) where larger genome and result files will be stored. You can change the location of the data folder by modifying the `build.sh` script.

To login to the system (user *apauser*), simply run the `run_apauser.sh` script. You are now running the Docker container with all required dependencies and software to run the `apa` example provided.

### Example

To directly run the example, skip to [Running the example](#running-the-example).

We provide an example Lexogen Quantseq Reverse sequencing run consisting of 6 experiments (3 HEK293, 3 TDP-43 KD) from the publication:

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />

In total 6 FASTQ files. The reads will be downloaded and processed in the Docker container by the `docker/example.sh` script. Reads for the example were preselected to match only the ones mapping to the chr22 of the hg38 assembly.

#### Data folder organization and structure

Data is organized inside **libraries**. Each library contains several experiments. Each experiment is represented by a FASTQ file. The sequencing (library) data is stored in the `apa.path.data_folder`. The location of the data folder can be changed by editing `apa/config/config.txt` and adding a line with `apa.path.data_folder="/path/to/data_folder"`.

```
# structure the data folder on apa

->data_folder [folder]
  ->example [folder] # folder of sequencing library with id example
    ->annotation.tab [file] # annotation TAB separated file (see below for structure)
    ->example.config # config file for the entire library (see below for structure)
    ->e1 [folder] # experiment 1
      ->example_e1.fastq.gz [file] # FASTQ file of experiment 1
      ->m1 [folder] # mapping 1, usually, we only map each experiment once, however several mappings (diverse parameters) can be addded (m1, m2, ...)
    ->e2 [folder] # experiment 2
      ->example_e2.fastq.gz [file] # FASTQ file of experiment 2
  ->example2 [folder] # folder of sequencing library with id example2
  ->example3 [folder] # folder of sequencing library with id example3
  ...
```
##### Structure of annotation.tab

The `annotation.tab` file is present for each library and annotates the experiments within the library. It's a TAB separated file with the following structure:

```
# first line is always a comment
exp_id	species	map_to	method	condition	replicate
1	hg38chr22	hg38chr22	lexrev	HEK293	1
2	hg38chr22	hg38chr22	lexrev	HEK293	2
3	hg38chr22	hg38chr22	lexrev	HEK293	3
4	hg38chr22	hg38chr22	lexrev	KD	1
5	hg38chr22	hg38chr22	lexrev	KD	2
6	hg38chr22	hg38chr22	lexrev	KD	3
```

Each line represents an experiment.

```
exp_id = experiment id (integer number)
species = species, descriptive (human, mouse, ms, hg, etc.)
map_to = id of genome (linked to pybio package) for read mapping
method = id of method
condition = descriptive
replicate = descriptive
```

##### Structure of example.config

The minimum structure of the library config file is as follows:

```
name:HEK-239
notes:Library with 3'-end Lexogen Quantrev sequencing of HEK-293 cell lines
method:lexrev
genome:hg38chr22
map_to:hg38chr22
seq_type:single
tags:
status:
columns:[['Condition', 'condition'], ['Replicate', 'replicate']]
columns_display:[['Condition', 'condition'], ['Replicate', 'replicate']]
```

#### Running the example

To run the provided example, start `~/apa/docker/example.sh` inside the Docker container. This will download the chr22 of the hg38 genome assembly, download and map the 6 example experiments (3 HEK293 and 3 TDP-43 KD) to the hg38 genome (only chromosome 22). It will build a polyA database from the aligned reads, estimate read counts at the identified polyA sites and also provide gene expression.

## Installation and tryout as standalone

The best way to install `apa`, `pybio` and other dependencies on your own server is simply to follow the [Dockerfile](Dockerfile). In case of problems, open an Issue on this repository page.

### Dependencies

For a full list of dependencies, please refer to the [Dockerfile](Dockerfile).

* [pybio](https://github.com/grexor/pybio), install by following instructions on the GitHub page
* [R](https://www.r-project.org), install following instructions, recommended latest release or at least >4.0.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), `BiocManager::install("edgeR")`
* [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html), `BiocManager::install("DEXSeq")`
* [samtools](http://www.htslib.org/), `RUN pip3 install HTSeq`
* [regex](https://pypi.org/project/regex), `pip3 install regex`

## Authors

[apa](https://github.com/grexor/apa) is developed and supported by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2009 when Tomaž Curk and Gregor Rot wrote the first prototype of [apa](https://github.com/grexor/apa). In 2013, Gregor Rot refactored and further developed the code, also establishing [expressRNA](http://expressRNA.org), a web application for exploring results of alternative polyadenylation analysis.

## Citing apa

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/apa/issues) to report issues and leave suggestions.
