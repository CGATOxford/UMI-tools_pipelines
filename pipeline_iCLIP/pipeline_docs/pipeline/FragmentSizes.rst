Fragment sizes
===============

Fragment sizes are calculated in two different ways. If the fragment is shorter than the length of the read, minus the barcode and UMI, then it is given by the length of the trimmed read, remembering to substract the length of any intron in the read. If the read is not trimmed beyond the barcode and UMI, then if is larger than the read size and is given by the insert size (minus any intron in the first read in pair).

Raw mapped reads on the genome
-------------------------------



.. report:: Sample_QC.MappedFragLength
   :render: r-ggplot
   :transform: melt
   :groupby: all
   :statement: aes(x=track, y=value/sum(value), fill=variable) + geom_bar(position="fill", stat="identity") + ylab("Fraction of reads") + scale_fill_discrete(name="Length bin (bp)") + coord_flip() + theme_bw()

   Distribution of fragment sizes of the raw reads mapped on the genome


Deduped fragments
-------------------


.. report:: Sample_QC.DedupedFragLengths
   :render: r-ggplot
   :transform: melt
   :groupby: all
   :statement: aes(x=track, y=value/sum(value), fill=variable) + geom_bar(position="fill", stat="identity") + ylab("Fraction of reads") + scale_fill_discrete(name="Length bin (bp)") + coord_flip() + theme_bw()

   Distribution of fragment sizes after deduplication



Duplication rates by size
--------------------------

.. report:: Sample_QC.LengthDedupedRatios
   :render: interleaved-bar-plot
   

   Percent of reads unique after deduplication for each size category


