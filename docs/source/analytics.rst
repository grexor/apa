**********************************
Analytics (methodology of data analysis)
**********************************

.. _3protocols:

3' end sequencing protocols
===========================

apa-db.org supports several 3' end sequencing protocols. After read pre-processing and alignment to the reference, the main difference between
these protocols is in determining the cleavage site loci and in further filtering steps (see below).

======== ============================= =============
protocol description                   site loci
======== ============================= =============
pA-seq   Wang et al., unpublished      last nucleotide of alignment
lex_fwd  Lexogen 3' forward            last nucleotide of alignment
lex_rev  Lexogen 3' reverse            first nucleotide of alignment
======== ============================= =============

.. _r_bedgraph_method:

Mapping of reads to the reference genome
========================================

Reads mapping (alignment) to the appropriate reference genome (hg19, mm10, etc.) is done with `STAR <https://github.com/alexdobin/STAR/releases>`_.
No special preprocessing other than quality control is performed. `STAR <https://github.com/alexdobin/STAR/releases>`_ is run allowing 5' and 3' soft clipping:
the potential poly-A tail (poly-T in some protocols) is soft-clipped from the read, allowing a more accurate identification of cleavage sites compared
to pre-processing and removing A/T rich 3'/5' ends of reads prior to mapping.

.. figure:: figures/clipping_analysis.png
  :width: 900px
  :figwidth: 900px
  :align: center

  Clipping analysis of aligned reads. Green line shows percentage of aligned nucleotides at specific position, blue line clipping from 5' end of reads and red line clipping from 3' end of reads.


R (raw, unfiltered) bedGraph
============================

The raw :ref:`bedGraph <r_bedgraph_format>` sites file, constructed by providing one site loci per alignment. This information is extracted from the bam file.
The sites loci are determined based on :ref:`protocol <3protocols>` (5' or 3' end of alignment). Files are stored in:

.. code-block:: bash

 ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}.T.bg

T (tail, filtered) bedGraph
===========================
The tail :ref:`bedGraph <r_bedgraph_format>` sites file. The filtering depends on the protocol (see below).

pA-seq (Wang et al.)
#################################

Lexogen 3' forward
######################

The assignment of the cleavage site follows this procedure:

#. set the site loci (cs_loci) as the last (3' end) nucleotide of the alignment
#. check 18nt upstream of cs_loci: IF #A>10, skip alignment

Since we allow soft clipping and do not pre-process the reads, we need to check the ending of the alignment for internal priming.

Lexogen 3' reverse
######################

To-do docs

Local poly-A atlas (database)
===========================

Before computing expression files, we define local groups of experiments (usually every library has it's own poly-A atlas). The atlas definition file is stored in:

.. code-block:: bash

  ${polya_folder}/${lib_id}.config

This config file contains the experiment identifiers, e.g.:

.. code-block:: bash

  elib_e1
  elib_e2
  elib_e3
  elib_e4

We group together T bedGraph files from the defined experiments (in this example e1, e2, e3 and e4) and:

#. create sorted loci list (reverse order = most expressed loci at the top)
#. select first loci in the list and remove all surrounding loci in region [-125, 125]
#. store selected loci in the atlas and remove it from the list
#. repeat previous two steps until loci list is empty

E (expression) bedGraph
===========================

For a given poly-A atlas and experiment R file, compute:

#. sum up experiment R file in region [-100, 25] for each poly-A atlas loci
