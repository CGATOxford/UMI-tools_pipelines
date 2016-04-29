Distribution of UMIs
---------------------


.. report:: Sample_QC.UMI_stats
   :render: r-ggplot
   :layout: column-4
   :display: png,png,50
   :statement: aes(x=sample,y=freq) + geom_violin() + scale_y_log10() + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=1/(4^5), lty=2)

   Distribution of UMIs


Reads Per Sample
-----------------


.. report:: Sample_QC.ReadsPerSample
   :render: r-ggplot
   :layout: column-4
   :display: png,png,50
   :statement: aes(y=total, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x,...) format(x,...,big.mark=",", scientific= F, trim = T)) + ylab("Reads")

   Number of reads for each barcode


.. report:: Sample_QC.PercentDemuxed
   :render: r-ggplot
   :layout: column-4
   :display: png,png,40
   :statement: aes(y=demuxed * 100, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x)) + ylab("Percent Passed Filter")

   Reads succesfully demuxed and filtered


.. report:: Sample_QC.PercentMapped
   :render: r-ggplot
   :layout: column-4
   :display: png,png,40
   :statement: aes(y=mapped * 100, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x), limits = c(0,100)) + ylab("Percent reads mapped")

   Percent Reads mapped



.. report:: Sample_QC.PercentDeDuped
   :render: r-ggplot
   :layout: column-4
   :display: png,png,40
   :statement: aes(y=p_unique * 100, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x)) + ylab("Percent reads unique")

   Percent of Reads Unique



.. report:: Sample_QC.FinalReads
   :render: r-ggplot
   :statement: aes(y=reads_mapped, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x,...) format(x,...,big.mark=",", scientific= F, trim = T)) + ylab("Total unique mapped reads")

   Final total mapped reads

.. report:: Sample_QC.PercentSpliced
   :render: r-ggplot
   :statement: aes(Track, pspliced) + geom_bar(stat="identity") +  scale_y_continuous(labels = function(x) sprintf("%.0f%%",x*100)) + ylab("Percent reads spliced") + theme_bw() + theme(axis.text.x=element_text(angle=90))

   Percent of deduped reads spliced

.. report:: Sample_QC.ReadLengths
   :render: r-ggplot
   :transform: melt
   :statement: aes(x=track, y=value/sum(variable), fill=Slice) + geom_bar(position="fill", stat="identity") + ylab("Fraction of reads") + scale_fill_discrete(name="Length bin (bp)") + coord_flip() + theme_bw()

   Read length distribution of unmapped reads




Saturation Analysis
--------------------

.. report:: Sample_QC.AlignmentSaturation
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=subset, y=counts, color = factor, shape = factor) + geom_point() + geom_line() + facet_wrap(~replicate) + theme_bw() + theme(aspect.ratio = 1)

   Subsampling of alignments


.. report:: Sample_QC.AlignmentSaturation
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=counts/subset, y=counts, color = factor, shape = factor) + geom_point() + geom_line() + facet_wrap(~replicate) + theme_bw() + theme(aspect.ratio = 1)

   Tests for model assumptions


.. report:: Sample_QC.LibrarySize_Binom
   :render: r-ggplot
   :statement: aes(x=subset, y=alignments) + geom_point() + geom_line(aes(y=expected_unique)) + geom_hline(yintercept=rframe$library_size[1]) + theme_bw()
   :width: 200
   :layout: column-4
   

   curve fits for saturation using Binomal distribution



.. report:: Sample_QC.LibrarySize_mm
   :render: r-ggplot
   :statement: aes(x=subset, y=alignments) + geom_point() + geom_line(aes(y=expected_unique)) + geom_hline(yintercept=rframe$library_size[1]) + theme_bw()
   :width: 200
   :layout: column-4

   curve fits for saturation using reciprical fit


.. report:: Sample_QC.mm_fit_stats
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track, y=Library.Size) + geom_bar(stat="identity") + geom_bar(aes(y=Library.Size*Percent.Saturation/100), stat="identity", fill = "orange") + theme_bw() + theme(axis.text.x = element_text(angle=90))

   Library size estimates

Context Stats
---------------

.. report:: Sample_QC.ContextStats
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track, y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip() + ylab("fractraction alignmnets") + scale_fill_brewer(type="qual",palette="Paired")

   Mapping Context for dediped reads


.. report:: Sample_QC.ContextRepresentation
   :render: r-ggplot
   :statement: aes(category, log2(precent_alignments/percent_bases)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90,hjust=1)) + ylab("log2 enrichment")
   :layout: column-3
   :groupby: track

   Enrichments of contexted over expectation


Splicing Index
------------------

Splicing index is:

.. math:: SI = \log_2\frac{2\sum N^{Exon-Exon}}{\sum N^{Exon-Intron} + N^{Intron-Exon}}

only consitutative exons are used, and only reads that map exactly to both sides of the junction are counted.

.. report:: Sample_QC.SplicingIndex
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track, y=log2(SI)) + geom_bar(stat="identity") + coord_flip() + xlab("Track") + ylab("Splicing index")

   Splicing index for each track
