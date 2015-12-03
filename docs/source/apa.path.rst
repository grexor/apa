.. role:: green
.. raw:: html

  <style>
  .green {
    color:green;
  }
  </style>

**********************************
Path module (``apa.path``)
**********************************

The root_folder is the location of the installed apa platform and is determined at import time. You can then override the defaults and change the names of several sub-folders.

================ =========== ===========
variable         default     description
================ =========== ===========
data_folder      data.apa    library, experiment, mapping and bedGraph files
polya_folder     data.polya  poly-A atlas (database) files
comps_folder     data.comps  comparisons for searching of APA gene
iCLIP_folder     data.iCLIP  iCLIP data used for RNA-maps (in bedGraph format)
================ =========== ===========

.. autofunction:: apa.path.t_filename
.. autofunction:: apa.path.r_filename
