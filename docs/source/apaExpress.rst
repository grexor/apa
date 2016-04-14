**********************************
**apaExpress**
**********************************
In this section we start by shortly describing the use of the web application and continue with a more detailed description of apaExpress analysis concepts and methodology.

Using the web application
-------------------------
The main menu on the top right in the browser is your basic navigation point in apaExpress.

.. figure:: figures/help_menu.png
  :align: center

#. Search allows you to find comparisons and experiments in apaExpress. Use it to browse all sorts of results and interact with them.
#. Help opens this documentation.
#. About provides some compact information about apaExpress research platform, it's use and contributors.
#. Sign-in provides the possibility to login with your existing Google account for accessing un-published data on our servers at "apaexpress.org". We do not store any information or get access to any of your Google account details.

Search menu
===========
This is the main interface to find and interact with apaExpress analysis results. Currently, the main search is for :ref:`comparisons <comps>`.

.. _comps:

Comparisons
===========
The comparison is performed between a set of test experiments against a set of control experiments, with the goal of identifying alternatively polyadenylated genes.

.. figure:: figures/help_comps.png
  :align: center

#. Name of comparison
#. Short description of comparison
#. :ref:`RNA-map <rnamap>` of proximal polyA sites and splice-site 1 (below)
#. :ref:`RNA-map <rnamap>` of distal polyA sites and splice-site 2 (below)
#. Experimental setup: information about the included experiments (in test and control sets)
#. RNAmotifs: results of RNAmotifs analysis in 3 regions up/down of regulated polyA sites
#. GO analysis: GO-term enrichment analysis of regulated target genes

We separately consider :ref:`3 categories of APA <apacat>`: same-exon, combined-exon and skipped-exon.

.. _rnamap:

RNA-maps
--------

.. _apacat:

APA categories
--------------
Once we determine the 2 polyA sites of interest (proximal, distal) for each gene, we divide the analysis in 3 categories based on the polyA site pair position:

.. figure:: figures/apa_categories.png
  :align: center

#. same-exon: the 2 polyA sites are in the same exon
#. composite-exon: the proximal site is annotated to the intron
#. skipped-exon: the proximal site is in another exon compared to the distal site

We also show the splice-site1 and splice-site2 positions relative to the proximal and distal sites.

The 4 loci per gene (1=proximal, 2=distal, 3=splice-site1, 4=splice-site2) are used to draw RNA-maps.
