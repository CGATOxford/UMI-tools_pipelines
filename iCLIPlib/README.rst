Tools for dealing with iCLIP data
==================================

This package contains a selection of programming tools for dealing with iCLIP data
and a bunch of scripts based on them. 

There are several pipelines available for dealing with iCLIP data - turning reads into
significant bases, but often we any iCLIP project will require a number of non-standard 
analysis. Here we have a set of tools for doing ad-hoc analysis.

In general profiles are represented as pandas Series where the index represents 
genomic bases and the values represents tag counts. 

The principle functions that produce these are:
    * count_intervals
    * count_transcript

Useful tools for analysis are:
    * randomizeSites - randomise the sites in a profile
    * rand_apply - randomize a profile a number of times and apply a function to it
    * spread - take a profile and extend each tag some bases in each direction
    
Also useful is 
   * TranscriptCoordInterconverter - as class for converting genomic coordinates to transcript ones

In addition to this are implementations for a number of published algorythms:
   * pentamer_enrichment - for looking for enriched kmers compared to randomised profiles
   * meta_gene - calculate a meta gene profile 
   * get_crosslink_fdr_by_randomisation - find significantly crosslinked bases by comparison to randomised profiles
     
and a number of scripts that implement example analyses. 

The requirements are CGAT, pysam, numpy and pandas.
The pipeline has numerous other requirements. 

.. warning::
  This code is mostly undocumented at the moment. I'm also in the middle of refactoring the code.  It is untested outside of its developement environment. 



  

