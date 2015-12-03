.. role:: green
.. raw:: html

  <style>
  .green {
    color:green;
  }
  </style>

**********************************
Analytics (methodology of data analysis)
**********************************

.. _3protocols:

3' end sequencing protocols
===========================

apa-db.org supports several 3' end sequencing protocols. After read pre-processing and alignment to the reference, the main difference between
these protocols is in determining the cleavage site loci and in further filtering steps (see below).

======== ============================= =============
protocol description                   cleavage loci
======== ============================= =============
pA-seq   Wang et al., unpublished      last nucleotide of alignment
lex_fwd  Lexogen forward strand 3' end last nucleotide of alignment
lex_rev  Lexogen reverse strand 3' end first nucleotide of alignment
======== ============================= =============

.. _r_bedgraph_method:

Mapping of reads to the reference genome
========================================

Reads mapping (alignment) to the appropriate reference genome (hg19, mm10, etc.) is done with `STAR <https://github.com/alexdobin/STAR/releases>`_.
No special preprocessing other than quality control is performed. `STAR <https://github.com/alexdobin/STAR/releases>`_ is run allowing 5' and 3' soft clipping:
the potential poly-A tail (poly-T in some protocols) is soft-clipped from the read, allowing a more accurate identification of cleavage sites compared
to pre-processing and removing A/T rich 3'/5' ends of reads prior to mapping.

.. figure:: bamclip.png
  :width: 900px
  :figwidth: 900px
  :align: center

  Clipping analysis of aligned reads. Green line shows percentage of aligned nucleotides at specific position, blue line clipping from 5' end of reads and red line clipping from 3' end of reads.


R (raw, unfiltered) bedGraph
============================

The raw :ref:`bedGraph <r_bedgraph_format>` cleavage site files. Not filtered. One sincle nucleotide per read alignment extracted from the bam files.
The cleavage site position is determined based on the :ref:`protocol <3protocols>` (5' or 3' end of alignment). Files are stored in:

.. code-block:: bash

 ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}.T.bg

T (tail, filtered) bedGraph
===========================
The tail :ref:`bedGraph <r_bedgraph_format>` cleavage site files. The filtering depends on the protocol.

pA-seq (Wang et al.)
#################################

Lexogen forward 3' end
######################

some info

Lexogen reverse 3' end
######################

Some other info
