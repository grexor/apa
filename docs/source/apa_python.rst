============================
**apa** (python) software manual
============================

-----------
Quick start
-----------

.. raw:: html

  <b>apa</b> is the python module behind all computational analysis and results. <b>apa</b> implements several steps of the analysis:


  <ul>
  <li> organizes the experimental data into libraries and manages the experiment annotation
  <li> processes aligned reads and creates polyA atlases together with polyA expression
  <li> computes comparisons: identifying alternatively polyadenylated genes (APA, this is where the name comes from)
  <li> integrates RNA-protein binding data and draws RNA-maps
  </ul>

  In order to be able to compute all the above, <b>apa</b> also integrates several external software packages:

  <ul>
  <li> genomic annotation and low level processing: <a href="https://github.com/grexor/pybio" target=_pybio>pybio</a>
  <li> motif analysis: <a href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r20" target=_rnamotifs>RNAmotifs</a>
  <li> GO ontology enrichment analysis: <a href='http://orange-bioinformatics.readthedocs.io/en/latest' target=_orange_bio>Orange Bioinformatics Add-on</a>
  <li> differential gene expression: <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_edgeR">edgeR</a>
  <li> read alignment: <a href="https://github.com/alexdobin/STAR" target="_star">STAR</a>
  </ul>

The minimal set of external software to be installed for **apa** to operate is STAR and pybio. Go analysis requires Orange and it's Bioinformatics Add-on.

Minimal set of dependencies
---------------------------

Short instructions to install STAR and pybio:

.. code-block:: bash

  # STAR short-read aligner
  git clone https://github.com/alexdobin/STAR.git  # clone GitHub STAR repository
  export PATH=$PATH:STAR/bin/Linux_x86_64          # add the correct STAR binary to your PATH (in this example Linux_x86_64)

  # pybio genomic analysis
  git clone https://github.com/grexor/pybio.git pybio   # clone pybio repository
  export PATH=$PATH:pybio/bin                           # add pybio/bin to PATH
  # + download any genomes that you will be using with the provided .sh script in the pybio/genomes folder

Prepare your library annotation
-------------------------------
For this example we will call our library **elib** (unique library identifier). Let's assume it contains 4 experiments (with id 1-4).
Each experiment is represented by one FASTQ file. The libraries are stored in the ${data_folder} (apa/data by default):

.. code-block:: bash

  ${data_folder}/elib/annotation.tab        # annotation file describing the library experiments
  ${data_folder}/elib/e1/elib_e1.fastq.gz   # FASTQ for experiment 1
  ${data_folder}/elib/e2/elib_e2.fastq.gz   # FASTQ for experiment 2
  ${data_folder}/elib/e3/elib_e3.fastq.gz   # FASTQ for experiment 3
  ${data_folder}/elib/e4/elib_e4.fastq.gz   # FASTQ for experiment 4

Example annotation.tab (TAB delimited file):

====== ======= === ====== ======= ======= ======
exp_id method  rep tissue cond    species map_to
====== ======= === ====== ======= ======= ======
1      lexfwd 1   HeLa   control Hs      hg19
2      lexfwd 2   HeLa   control Hs      hg19
3      lexfwd 1   HeLa   test    Hs      hg19
4      lexfwd 2   HeLa   test    Hs      hg19
====== ======= === ====== ======= ======= ======

It contains 4 experiments, 2 test, 2x control, on HeLa cells, sequenced with Lexogen forward 3' end method.


Define comparison
-----------------

We need to define which test experiments are going to be compared to which controls. We define an elib.config file
in the ${comps_folder} (by default apa/data_comps):

.. code-block:: bash

  ${comps_folder}/elib/elib.config   # elib comparison configutation file

In elib.config (TAB delimited file), we assign our experiments to control or test:

.. code-block:: bash

  id	experiments	name
  c1	elib_e1	    control_1
  c2	elib_e2	    control_2
  t1	elib_e3	    test_1
  t2	elib_e4	    test_2

At this point, we have a library with experiments in the ${data_folder}/elib and one comparison in the ${comps_folder}/elib.

Process data and compute comparison
-----------------------------------

First we map the reads to the reference genome (hg19):

.. code-block:: bash

  apa.map -comps_id elib

Then we generate bedGraph files from the aligned reads (depending on the protocol, the first mapped nucleotide or the last mapped nucleotide is considered
and the strand is reversed or kept):

.. code-block:: bash

  apa.bed.multi -lib_id elib

Then we create the polyA atlas from all the hg19 experiments:

.. code-block:: bash

  apa.polya -poly_id hg19

and finally we take the final atlas of polyA sites to compute expression in each of the experiments (results are again stored in bedGraph files):

.. code-block:: bash

  apa.bed.multi -lib_id -type expression

Now we have the atlas of polyA sites in our data and also the expression of the sites in experiments 1 to 4. We can analyze the data and compute the comparison:

.. code-block:: bash

  apa.comps -comps_id elib

Comparison results are stored in the ${comps_folder}/elib.

-------
Methods
-------

.. 3_protocols:

3' end sequencing protocols
-------------------------------

apaExpress supports several 3' end sequencing protocols. Basically we could divide the procedure of processing reads from different sequencing protocols in two
categories: (1) the polyA site (cleavage site) position is at the first aligned nucleotide of each aligned read, or (2) at the last nucleotide. Depending on the
protocol, the strand could be either left the same (as the strand of the aligned read) or reversed.

=========================================== ================================ ===================
Protocol                                    Cleavage site                    Strand
=========================================== ================================ ===================
fwd (e.g. Lexogen 3' forward, pA-seq)       3'-end nucleotide of alignment   leave or reverse
rev (e.g. Lexogen 3' reverse, PolyA-seq)    5'-end nucleotide of alignment   leave or reverse
=========================================== ================================ ===================

.. read_alignment:

Mapping of reads to the reference genome
-------------------------------

.. raw:: html

  Reads mapping (alignment) to the reference genome (hg19, mm10, etc.) is done with <a href=https://github.com/alexdobin/STAR/releases target="_star">STAR</a>.
  No special preprocessing other than quality control is performed. STAR is run allowing 5' and 3' soft clipping: the potential poly-A tail (poly-T in some protocols)
  is soft-clipped from the read, allowing a more accurate identification of cleavage sites compared to pre-processing and removing A/T rich 3'/5' ends of reads
  prior to mapping.

.. figure:: figures/star_clipping.png
  :width: 600px
  :align: center

  The above figure shows the amount of clipping at 5' and 3' end of aligned reads.

Custom polyA atlas
-------------------------------

Instead of constructing the polyA atlas from all experiments available (species wise), it's possible to define a custom set of experiments for the polyA atlas.
The custom set is a single file, e.g. an atlas with name "myatlas" would define experiments in the file ${polya_folder}/myatlas.config:

.. code-block:: bash

  ${polya_folder}/myatlas.config

This config file contains the experiment identifier per line, e.g.:

.. code-block:: bash

  elib_e1
  elib_e2
  elib_e3
  elib_e3

------------
File formats
------------

Description of various file formats with their structure that apa platform supports.

bedGraph files
--------------
The format is classic bedGraph: chr, strand, [pos_start, pos_end). An example bedGraph:

.. code-block:: python
  :linenos:

  track type=bedGraph name="e1" description="HeLa cells" db=hg19 color="120,101,172" priority="20" maxHeightPixels="100:50:0" altColor="200,120,59" visibility="full"
  chr1  1200  1230  100
  chr1  2000  2100  -30

Line 1 contains the track information (can be omitted). Line 2 defines a region on chr1 from loci 1200 to 1229 with value 100 (positive strand).
Line 2 defines another region on the same chromosome from 2000 to 2099, this time on the negative strand (value -30).

For a detailed description of bedGraph files, see the `UCSC description <http://genome.ucsc.edu/goldenpath/help/bedgraph.html>`_.

--------------------------------------
Annotation module (``apa.annotation``)
--------------------------------------

.. autofunction:: apa.annotation.read

------------------------
Bed module (``apa.bed``)
------------------------

.. autofunction:: apa.bed.write_bed
.. autofunction:: apa.bed.bed_raw
.. autofunction:: apa.bed.bed_expression

--------------------------
Path module (``apa.path``)
--------------------------

The root_folder is the location of the installed apa platform and is determined at import time. You can then override the defaults and change the names of several sub-folders.

================ =========== ===========
variable         default     description
================ =========== ===========
data_folder      data.apa    library, experiment, mapping and bedGraph files
polya_folder     data.polya  poly-A atlas (database) files
comps_folder     data.comps  comparisons for searching of APA gene
iCLIP_folder     data.iCLIP  iCLIP data used for RNA-maps (in bedGraph format)
================ =========== ===========
