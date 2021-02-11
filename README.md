# apa: alternative polyadenylation (APA) analysis

* [About](#about)
* [Installation and tryout with Docker](#installation-and-tryout-with-docker)
* [Installation and tryout as standalone](#installation-and-tryout-as-standalone)
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

## Installation and tryout with Docker

Clone this repository (`git clone https://github.com/grexor/apa.git`) and run `build.sh`.

This will build a Docker image with *apa*, *pybio* and all other dependencies installed. It will also create a "data" folder on your local drive (inside the same folder where build.sh is) where larger genome and result files will be stored. You can change the location of the data folder by modifying the `build.sh` script.

To login to the system (user *apauser*), simply run the `run_apauser.sh` script. To run the provided example, run `apa/docker/example.sh` inside the Docker container. This will map the 6 example experiments to the hg38 genome (only chromosome 22).

## Installation and tryout as standalone

### Clone the GitHub repository

Clone this repository and add the containing folder to PYTHONPATH:

```
git clone https://github.com/grexor/apa.git
```

If, for example, you installed `apa` to `/home/user/apa`, you would add this command to the `.profile` file in your home folder:

```
export PYTHONPATH=$PYTHONPATH:/home/user
export PATH=$PATH:/home/user/apa/bin
```

### Dependencies

* [pybio](https://github.com/grexor/pybio), install by following instructions on the GitHub page
* [R](https://www.r-project.org), install following instructions, recommended latest release or at least >4.0.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), `BiocManager::install("edgeR")`
* [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html), `BiocManager::install("DEXSeq")`
* [samtools](http://www.htslib.org/), `RUN pip3 install HTSeq`
* [regex](https://pypi.org/project/regex), `pip3 install regex`

For a full list of dependencies, please refer to the [Dockerfile](Dockerfile).

## Documentation

Here we provide basic `apa` usage examples.

### Download and prepare example library

For this example we will call our library `example` (unique library identifier). Let's assume it contains 6 experiments (with ids e1-e6). Each experiment is represented by one FASTQ file. The libraries are stored in the `apa.path.data_folder` (`apa/data.apa` by default):

```
apa/data.apa/example/example.config          # config file where you define species, reference genome and sequencing protocol
apa/data.apa/example/annotation.tab             # annotation file describing the experiments of data library 20201104_1
apa/data.apa/example/e1/example_e1.fastq.gz     # FASTQ for experiment 1
apa/data.apa/example/e2/example_e2.fastq.gz     # FASTQ for experiment 2
apa/data.apa/example/e3/example_e3.fastq.gz     # FASTQ for experiment 3
apa/data.apa/example/e3/example_e4.fastq.gz     # FASTQ for experiment 4
apa/data.apa/example/e4/example_e5.fastq.gz     # FASTQ for experiment 5
apa/data.apa/example/e4/example_e6.fastq.gz     # FASTQ for experiment 6
```

You can download the files from the expressRNA.org server by running:

```
cd apa/data.apa/example
for exp_id in 1 2 3 4 5 6
do
  mkdir e${exp_id}
  wget https://expressrna.org/share/data/example/e${exp_id}/example_e${exp_id}.fastq.gz -O e${exp_id}/example_e${exp_id}.fastq.gz
done
```

### Process example library

To map (align) the example library 4 fastq files to the reference genome (hg38), just run:

```
apa.map.lib -lib_id example
```

Next, run these commands to process the library:

```
apa.bed.multi -lib_id example           # creates bedGraph files from the mapped reads
mkdir apa/data.polya                    # create folder to contain polyA atlas
apa.polya.makeconfig -lib_id example    # create config file to include all experiments in the polyA atlas
apa.polya -poly_id example              # create polyA atlas from all experiments in the library
```

Next, we use the created polyA atlas to compute "expression" (counts) for the polyA sites for each experiment:

```
apa.bed.multi -lib_id example -type expression -poly_id example -upstream 10 -downstream 10
```

### Prepare and process example control vs test comparison

Coming soon

## Authors

[apa](https://github.com/grexor/apa) is developed and supported by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2009 when Tomaž Curk and Gregor Rot wrote the first prototype of [apa](https://github.com/grexor/apa). In 2013, Gregor Rot refactored and further developed the code, also establishing [expressRNA](http://expressRNA.org), a web application for exploring results of alternative polyadenylation analysis.

## Citing apa

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/apa/issues) to report issues and leave suggestions.
