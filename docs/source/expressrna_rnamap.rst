.. _rnamap:

RNA-maps
--------

.. raw:: html

  RNA-maps are scatterplots of RNA-protein binding information around regulated polyA sites and relevant splice-sites in their surrounding. The maps show
  where and how much the protein of interest binds and also in which case (<font color=blue>repressed sites = blue</font>, <font color=red>enhanced sites = red</font>).

.. figure:: figures/rna_maps.png
  :align: center
  :width: 750px
  :figwidth: 750px

#. Name of CLIP data used to draw the map
#. Type of proximal, distal, splice-site1 and splice-site2 that the map displays (:ref:`see APA categories<apacat>`)
#. RNA-map of proximal polyA site. Blue is percentage of regulated genes (targets) and in red percentage of enhanced genes bound by the protein of interest at specific location in the **range [-200..200]** with polyA site at the center [0]
#. Mouse over displays more exact positional information, e. g. at position -69 (relative to polyA site), 3.5% of repressed proximal polyA sites are bound by the protein of interest
