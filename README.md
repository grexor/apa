# apa: alternative polyadenylation (APA) analysis

* [About](#about)
* [Installation](#installation)
* [Authors](#authors)
* [Reporting problems](#reporting-problems)

## About

apa is a Python framework for processing and analysing 3'-end targeted sequence data to study alternative polyadenylation. [apa](https://github.com/grexor/apa) interconnects [pybio](https://github.com/grexor/pybio) (basic handling of annotated genomes), [RNAmotifs2](https://github.com/grexor/rnamotifs2) (analysis of regulatory motif clusters) and other open-source software (DEXSeq, STAR short-read aligner).

The inclusive nature of the framework, together with novel integrative solutions (differential polyA site usage and RNA-protein binding via RNA-maps, cluster motif analysis), results in the following computational capabilities:

+ management of diverse high-throughput sequencing datasets (pre-processing, alignment, annotation),
+ polyA site database (atlas) construction and comparison to existing polyA resources,
+ identification of genes that undergo alternative polyadenylation (DEXSeq),
+ identification of motifs influencing polyA site choice (RNAmotifs2),
+ identification of motifs influencing alternative splicing (DEXSeq and RNAmotifs2),
+ integration with iCLIP (RNA-protein binding) and computing RNA-maps,
+ and other.

## Installation

A few steps of how to download and setup `apa`.

### Clone the GitHub repository

For now the most direct way of installing pybio is to clone the repository and add the containing folder to PYTHONPATH:

```
git clone https://github.com/grexor/apa.git
```

If, for example, you installed pybio to /home/user, you would add this command to the .profile file in the same folder:

```
export PYTHONPATH=$PYTHONPATH:/home/user
```

### Dependencies

There are a few software tools pybio depends on:

* [pybio](https://github.com/grexor/pybio), install by following instructions on the GitHub page
* [R](https://www.r-project.org), install following instructions, recommended latest release or at least >4.0.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), `BiocManager::install("edgeR")`
* [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html), `BiocManager::install("DEXSeq")`

## Authors

[apa](https://github.com/grexor/apa) is developed and supported by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2009 when Tomaž Curk and Gregor Rot wrote the first prototype of [apa](https://github.com/grexor/apa). In 2013, Gregor Rot refactored and further developed the code, also establishing [expressRNA](http://expressRNA.org), a web application for exploring results of alternative polyadenylation analysis.

## Citing apa

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/apa/issues) to report issues and leave suggestions.
