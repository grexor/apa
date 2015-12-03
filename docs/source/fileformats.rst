.. role:: green
.. raw:: html

  <style>
  .green {
    color:green;
  }
  </style>

**********************************
File formats
**********************************

Description of various file formats with their structure that apa platform supports or generates.

.. _r_bedgraph_format:

bedGraph files
==============
The format is classic bedGraph: chr, strand, [pos_start, pos_end). An example bedGraph:

.. code-block:: python
  :linenos:

  track type=bedGraph name="e1" description="HeLa cells" db=hg19 color="120,101,172" priority="20" maxHeightPixels="100:50:0" altColor="200,120,59" visibility="full"
  chr1  1200  1230  100
  chr1  2000  2100  -30

Line 1 contains the track information (can be omitted). Line 2 defines a region on chr1 from loci 1200 to 1229 with value 100 (positive strand).
Line 2 defines another region on the same chromosome from 2000 to 2099, this time on the negative strand (value -30).

For a detailed description of bedGraph files, see the `UCSC description <http://genome.ucsc.edu/goldenpath/help/bedgraph.html>`_.
