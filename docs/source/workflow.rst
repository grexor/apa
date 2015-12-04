**********************************
Quick start
**********************************

A quick start example to show you how to use the platform to process your data.

Pre-requisites
-------------------------------

This list of software needs to be installed and in your path to run the analysis:

.. code-block:: bash

  # STAR short-read aligner
  git clone https://github.com/alexdobin/STAR.git  # clone GitHub STAR repository
  export PATH=$PATH:STAR/bin/Linux_x86_64          # add the correct STAR binary to your PATH (in this example Linux_x86_64)

  # pybio genomic analysis
  git clone https://github.com/grexor/pybio.git pybio   # clone pybio repository
  export PATH=$PATH:pybio/bin                           # add pybio/bin to PATH

Download library and annotation
-------------------------------
Experiments are organized in libraries. For this example we will call our library example_lib, containing 4 experiments (with id 1, 2, 3 and 4).

To download the example library (elib) with 4 experiments (id 1,2,3 and 4), move to your ${data_folder} (by default data.apa) and run wget:

.. code-block:: bash

  wget http://www.apa-db.org/example . -R

This will download the experiment FASTQ files (one file per experiment) and the annotation file to your ${data_folder}:

.. code-block:: bash

  ${data_folder}/example_lib/annotation.tab # annotation file
  ${data_folder}/elib/e1/elib_e1.fastq.gz   # FASTQ experiment 1
  ${data_folder}/elib/e2/elib_e2.fastq.gz   # FASTQ experiment 2
  ${data_folder}/elib/e3/elib_e3.fastq.gz   # FASTQ experiment 3
  ${data_folder}/elib/e4/elib_e4.fastq.gz   # FASTQ experiment 4

Check annotation.tab:

====== ======= === ====== ======= ======= ======
exp_id method  rep tissue cond    species map_to
====== ======= === ====== ======= ======= ======
1      lex_fwd 1   HeLa   control Hs      hg19
2      lex_fwd 2   HeLa   control Hs      hg19
3      lex_fwd 1   HeLa   test    Hs      hg19
4      lex_fwd 2   HeLa   test    Hs      hg19
====== ======= === ====== ======= ======= ======

It contains 4 experiments, 2x test and 2x control, from HeLa cells.

Map reads to reference genome (hg19)
------------------------------------

Genomic data processing is done with pybio. Download and setup the hg19 reference genome and Ensemble annotation:

.. code-block:: bash

  cd pybio/genomes     # cd to genomes folder
  ./hg19.download.sh   # download hg19 assembly and Ensembl annotation, build STAR indices
