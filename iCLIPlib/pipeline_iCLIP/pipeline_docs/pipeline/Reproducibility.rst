Reproducibility
----------------

Reproducilbity measures the number of sites with at least n reads mapping to them in one replicate that have reads mapping to them in 1 or 2 of the other replicates as a fraction of the total number of sites with that depth in that replicate. 

.. report:: Sample_QC.Reproducibility
   :render: r-ggplot
   :groupby: all
   :statement: aes(level, reproducibility, color=Replicate) + geom_line() + geom_point() + facet_grid(slice ~ track) + coord_cartesian(xlim=c(0,5)) + theme_bw()
   :tf-label-level: 3

   Reproduciblity


The problem with the measure above (which is the one outlined in Sutomui et al) is that the largest rep will always have a lower reproducibility because all those extra locations can't possibly be replicated. Below I normalise the reproduciblility by the maximum possible level of reproduciblitity.

.. report:: Sample_QC.NormReproducibility
   :render: r-ggplot
   :groupby: all
   :statement: aes(level, reproducibility, color=Replicate) + geom_line() + geom_point() + facet_grid(slice ~ track) + coord_cartesian(xlim=c(0,5), ylim=c(0,1)) + theme_bw()
   :tf-label-level: 3

   Normalised Reproduciblity

The next plot shows how reproducible cross-linked bases are in the control samples rather than in other replicates of the sample cell line. 

.. report:: Sample_QC.ReproducibilityVsControl
   :render: r-ggplot
   :groupby: all
   :slices: 1,3
   :statement: aes(level,reproducibility,color=Replicate) + geom_line() + geom_point() + facet_grid(slice~track) + theme_bw() + coord_cartesian(xlim = c(0,25))
   :tf-label-level: 3

   Reproducibility vs. Controls


Given that there is some reproducibility between one replicate and others pulling down the same factors and also some between that same replicate and the negative controls, how much infomation is there in the sample that is due to the correct pull down. Assuming that infomation shared between a sample and a control will also be shared by another replicate of the same pull down, the ratio of replicating bases between A) a replicate and a the control and B) one replicate and another should be above one, and the excess should speak to how much extra, factor specific information there is. 


.. report:: Sample_QC.ReproducibilityReplicateVsControl
   :render: r-ggplot
   :groupby: all
   :statement: aes(depth,ratio,color=Replicate) + geom_line() + geom_point() + facet_grid(slice~track) + scale_y_log10() + coord_cartesian(xlim=c(0,10)) + theme_bw()
   :tf-label-level: 3
   :slices: 1

   Ratio of reproducibility in replicates of same factor to that in other factors.



The reproducibility can also be used to calculate a distance metric between samples. The jaccard index is the interection of two sets divided by the union. By applying this accross each pair of samples at the 1 level we can build a clustering of samples.

