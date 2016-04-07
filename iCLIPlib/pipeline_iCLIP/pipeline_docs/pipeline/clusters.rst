Cluster Calling
================

Tracks of the clusteres can be found `here <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgct_customText=https://www.cgat.org/downloads/N6Cduavf7p/iCLIP_fullrun2/export/clusters/UCSC.txt>`_ 

Numbers of clusters called:


.. report:: clusters.ClusterCounts
   :render: r-ggplot
   :statement: aes(paste(sample,replicate), count, fill = replicate == "reproducible") + geom_bar(stat="identity", position="dodge") + theme_bw() + theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1)) + xlab("") + scale_y_continuous(name = "Clusters", labels=function(x) format(x, scientific = F, big.mark=","))

   Numbers of significant clusters called in each replicate and those reproducilbe


.. report:: clusters.ClusterCounts
   :render: table
   
   Numbers of significant clusters called in each replicate and those reproducilbe
   

.. _clustercontext:

Cluster Contexts
-------------------

.. report:: clusters.ClusterContexts
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track, y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip()

   Genomic contexts of significant clusters
