.. _pipelinemotifs:

Motif Finding
==============

MEME motifs within clusters
-----------------------------

Reproducilbe clusters
++++++++++++++++++++++

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: track
   :tracks: r(reproducible)
   :force:

   Top 5 motifs found in clusters from each track


Indevidual Replicates
++++++++++++++++++++++

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: track
   :tracks: r(.+-.+-R[0-9]+)
   :force:

   Top 5 motifs found in clusters from each track


DREME motifs within clusters
----------------------------

Reproducible clusters
++++++++++++++++++++++

.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :tracks: r(reproducible)
   :large-html-class: sortable
   :groupby: track

   Over-represented motifs in reproducible clusters compared to randomized clusters


Indevidual Replicates
++++++++++++++++++++++

.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :tracks: r(.+-.+-R[0-9]+)
   :large-html-class: sortable
   :groupby: track

   Over-represented motifs in reproducible clusters compared to randomized clusters
