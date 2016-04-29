'''
find_significant_bases.py - template for CGAT scripts
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script calculates p values for cross linked positions in transcripts from
an iCLIP experiment. Input (from stdin) is a gtf file containing transcripts
and (specified to the script) a bamfile containing the mapped reads. The output
is a bedgraph file containing the (optionally corrected) pvalues for
significant cross-link clustering around each crosslinked base.

The model is based on several assumptions,

1)  The null probability of a cross-link at any position in a transcript is
    equal. Thus:

2)  The probability of getting x_i read counts at any given position is
    distributed binomially. :math:`X_i ~ Bin(n,1/p)` where n is the total
    number of counts mapping to a region and p is the number of bases in the
    region.

3)  The probability of getting x_win read counts within a window around a
    cross-linked base is distributed :math:`X_win ~ Bin(n-x_i, 1/p_win)` where
    p_win is p minus the window size.

4)  Thus the probability of getting x_i reads or more at position i and X_win
    or more reads around it is:
    
    :math:`P(X_i>=x_i)*P(X_win>=x_win)`

5)  All transcripts from a gene are equally expressed. Where a crosslink
    position is with in more than one transcript the probability for that base
    is. Thus for 2 transcripts:
    
    :math:`P(X>x) = sum (P(X>=x| transcript i)*P(transcript i))`
           for i in transcripts

    as all transcripts are equally expressed this becomes

    :math:`P(X>x) = mean(P(X>=x | transcript i))`
           for i in transcripts.

    i.e. overlapping bases are averaged. This is probably not good. If you
    don't want this please deal with overlapping transcripts before running
    the script.

.. Note:: Input GTF should be sorted such that all lines from the same
         transcript are consecutive, and all transcripts from the same gene are
         consecutive.

.. Warning:: Genes must not overlap. If they do the overlapping bases will
             appear more than once in the output.

.. TODO:: Allow normalisation by rnaseq input
Options
-------

-g, --groupby: How to divide the transcript into regions that will be tested
               together. The default is exons. Here all exons from a transcript
               are concatenated, and treated as a single region, and the
               introns likwise.

               Alternatively, the utrs can also be treated as seperate regions
               (--groupby=utrs). If this is the case the GTF must contain CDS
               entires. Or the whole genomic region from tss to tts can be
               treated as a single region.

-w, --window-size: The script uses windows around the cross-linked base to
               determine clustering. This parameter specifies the distance
               either side of the cross-linked base to use.

-p, --pipeout: Specifies that results should be written to the output as soon
               as generated. This will save memory, and allow downstream
               computations, but may produce a truncated output if there is an
               error and does not allow an FDR correction.

-f, --fdr:     Compute an BH FDR correction on the results.
               Implies not --pipeout.

-t, --dtype:   The numpy dtype to use for storing counts. The default is
               uint32. Smaller types will use less memory, but run the risk of
               integer overflow (detected).

-
Usage
-----

.. Example use case

Example::

   python find_significant_bases.py mybam.bam > outfile.bg < transcripts.gtf

Usage::

   python find_significant_bases.py [OPTIONS] BAMFILE < GTFFILE


Command line options
--------------------

'''

from scipy.stats import binom
import pandas as pd
import numpy as np
import sys
import CGAT.Experiment as E
import iCLIP
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import pysam
from statsmodels.stats.multitest import multipletests
import CGAT.Intervals as Intervals
import CGAT.Bed as Bed


def get_windows(pvalues, window_size, threshold):

    # intervals are close closed
    windows = [(pos-window_size, pos+window_size+1)
               for pos in pvalues.index.values]

    merged_windows = Intervals.combine(windows)
    windows_min_p = [pvalues.ix[float(start):float(end-1)].min()
                    for start, end in merged_windows]
    return zip(merged_windows, windows_min_p)


def windows2bed12(windows, contig, strand, name, score):
    '''Convert a list of intervals into a single bed12 entry '''

    windows = sorted(windows)

    entry = Bed.Bed()

    #if strand == "-":
    #    windows = [(y+1, x+1) for x, y in windows]
    #    windows = sorted(windows)
    #else:
    #    windows = sorted(windows)

    entry.start = int(windows[0][0])
    entry.end = int(windows[-1][1])

    entry.contig = contig
    
    blockCount = int(len(windows))
    blockSizes = ",".join([str(int(window[1]-window[0])) for window in windows])
    blockStarts = ",".join([str(int(window[0] - windows[0][0])) for window in windows])
    thickStart = int(entry.start)
    thickEnd = int(entry.end)
    itemRGB = "255,0,0"
    
    entry.fields = [name, score, strand, thickStart, thickEnd, itemRGB,
                    blockCount, blockSizes, blockStarts]

    assert entry.end - entry.start > 0, "Malformed Bed entry entry size less than zero"
    assert all([blockSize > 0 for blockSize in 
                map(int, blockSizes.split(","))]), \
        "Malformed Bed entry, at least one block size less than zero"
    assert all([entry.start + blockStart <= entry.end
                for blockStart in map(int,blockStarts.split(","))]), \
                    "Malformed Bed entry: block start after end of entry"

    return entry
    

def bases_to_windows(pvalues, gene, window_size, threshold):

    contig = gene[0][0].contig
    strand = gene[0][0].strand
    gene_id = gene[0][0].gene_id
    try:
        gene_pvals = pvalues[gene_id][contig][strand]
    except KeyError:
        E.info("No Significant CLIP sites in gene %s" %
               (gene[0][0].gene_id))
        return []

    gene_pvals = gene_pvals.sort_index()
    outlist = []

    for transcript in gene:

        coords_converter = iCLIP.TranscriptCoordInterconverter(transcript)

        # first exons
        exons = GTF.asRanges(transcript, "exon")
       
        # pandas indexing is inclusive, but our exon intervals are
        # half closed
        exons_pvals_list = [gene_pvals.ix[float(start):float(end-1)]
                            for start, end in exons]

        exons_pvals = pd.concat(exons_pvals_list)

        exons_pvals.index = coords_converter.genome2transcript(
            exons_pvals.index.values)

        windows = get_windows(exons_pvals, window_size, threshold)
        windows = [(coords_converter.transcript_interval2genome_intervals(
            window),p) for window,p in windows]

        # now for introns
        for intron in GTF.toIntronIntervals(transcript):
            intron_pvals = gene_pvals.ix[float(intron[0]):float(intron[1]-1)]
            intron_windows = get_windows(intron_pvals, window_size, threshold)
            intron_windows = [((max(intron[0], start), min(intron[1], end)),p)
                              for (start, end), p in intron_windows]
            windows.extend([([window],p) for window,p in intron_windows])
        
        try:
            outlist.extend(
                [windows2bed12(window, contig, transcript[0].strand,
                               "%s_%s" % (transcript[0].transcript_id, n),
                               score=p)
                 for n, (window,p) in enumerate(windows)])
        except:
            print [x for x in enumerate(windows)]
            raise
        assert not any([bed.end - bed.start < 1 for bed in outlist])

    return sorted(outlist, key=lambda x: x.start)

        
class DeferredOutput:
    ''' This class looks like a file like object, but stores up all
     the objects passed for output, and outputs them all when close
    is called. Optionally performs multiple testing correction '''

    def __init__(self, outfile_bases=None, outfile_windows=None,
                 correct=False, window_size=0, threshold=0.05):

        self.outfile_bases = outfile_bases
        self.outfile_windows = outfile_windows
        self.output = []
 
        self.correct = correct
        self.window_size = window_size
        self.threshold = threshold
        self.genes = []
        self.max_end=0

    def write(self, gene_results, gene):


        self.output.append(gene_results)
        if self.outfile_windows:
            self.genes.append(gene)

    def close(self):

        output = pd.concat(self.output)
        output = output.sort_index()
        E.debug("most 3' coordingate seen is %s" % (output.index.values[-1],))
        if self.correct:
            E.info("Correcting p-values using BH ...")
            corrected_pvals = multipletests(output, method="fdr_bh")
            output = pd.Series(corrected_pvals[1], index=output.index)

        E.info("Writing output")
        E.debug("output contains %i entries" % len(output))
 
        if self.outfile_windows:
            E.info("Writing windows")
            sig_windows = output[output < self.threshold]
            for gene in self.genes:
                windows = bases_to_windows(sig_windows, gene, self.window_size,
                                           self.threshold)
                for bed in windows:
                    self.outfile_windows.write(str(bed) + "\n")

        if self.outfile_bases:
            E.info("Writing bases")
           
            output = output.reset_index()
            output.drop("strand", axis=1, inplace=True)
            output.drop("gene_id", axis=1, inplace=True)
            output = output.groupby(["contig", "position"], as_index=False).min()
            output.position = output.position.astype("int64")
            output["end"] = output["position"] + 1
            output = output[["contig", "position", "end", 0]]

            output.to_csv(self.outfile_bases,
                          sep="\t",
                          header=False,
                          index=False)


class InstantOutput:
    ''' This class looks file a file like object but takes pandas Series
    objects and outputs them to a file handle it keeps open '''

    def __init__(self, outfile_windows=None, outfile_bases=None,
                 window_size=0, threshold=0.05, **kwargs):

        self.window_size = window_size
        self.threshold = threshold
        self.outfile_windows = outfile_windows
        self.outfile_bases = outfile_bases

    def write(self, gene_results, gene):

        if self.outfile_bases:
            gene_results = gene_results.sort_index()
            pvalues = gene_results[gene_results < self.threshold]
            windows = bases_to_windows(pvalues, gene,
                                       self.window_size, self.threshold)
            for bed in windows:
                self.outfile_windows.write(str(bed))

        if self.outfile_windows:
            gene_results = pd.DataFrame({"start": gene_results,
                                         "end": gene_results})
            gene_results.reset_index()
            gene_results.to_csv(self.outfile_bases,
                                sep="\t",
                                header=False,
                                index=False)

    def close(self):
        
        try:
            self.outfile_bases.close()
        except AttributeError:
            pass

        try:
            self.outfile_windows.close()
        except AttributeError:
            pass


def calculateProbabilities(counts, window_size, length, start=0):
    '''Calculates the probablity of observing the counted
    number of reads in windows of "window_size" around each
    cross-linked size.

    Currently assumes that coordinates are in transcript space.

    The length of the transcript must be provided because it
    is needed for calculating the prior and may be outside of the
    provided coordinates.

    Start allows, together with length, for only using part of
    counts.'''

    # limit to subset
    counts = counts[(counts.index.values >= start) &
                    (counts.index.values < (start + length))]
    total_counts = counts.sum()
   
    # probability is counts-2 because we want P(X>=x) which is
    # 1 - P(X<x-1). Thats -1. The other -1 comes from the fact
    # that we want the p that any base in the transcript has
    # X>=x, not just this specific one.

    single_base_ps = 1 - binom.cdf(counts-2, total_counts, 1.0/length)

    heights = np.zeros(len(counts))

    window_start = np.maximum(0, counts.index.values-window_size)
    window_end = np.minimum(start+length, counts.index.values + window_size)

    ps = (window_end - window_start + 0.0) / length
    ps = ps.astype("float64")

    for i, base in enumerate(counts.index.values):

        try:
            window = counts[window_start[i]:window_end[i]]
        except KeyError:
            print (window_start, window_end)
            print counts
  
        heights[i] = window.sum()
       
    heights = heights - counts.values
    window_ps = pd.Series(1 - binom.cdf(heights - 1, total_counts, ps),
                          index=counts.index)

    # correct for number of independent windows.
    return window_ps*single_base_ps


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    grouping_choices = ["exons",
                        "utrs",
                        "all"]
    parser.add_option("-g", "--grouping", dest="grouping", type="choice",
                      choices=grouping_choices,
                      help="How to group transcript regions choices are [%s]"
                            % ",".join(grouping_choices))
    parser.add_option("-p", "--pipeout", dest="pipeout", action="store_true",
                      help="Output continuously to the pipe rather than in a"
                           "chunk at the end")
    parser.add_option("-d", "--dtype", dest="dtype", type="string",
                      default="int32",
                      help="Numpy dtype for storing counts")
    parser.add_option("-w", "--window-size", dest="window_size",
                      type="int", default=15,
                      help="Size of window either size of crosslinked base to"
                           "consider")
    parser.add_option("-f", "--fdr", dest="fdr", action="store_true",
                      default=False,
                      help="perform BH fdr correction on p-values, implies not"
                           "--pipeout")
    parser.add_option("-o", "--output-windows", dest="output_windows",
                      action="store_true",
                      default=False,
                      help="Output consolidated windows isntead of bases")
    parser.add_option("-b", "--output-both", type="string", dest="output_both",
                      default=None,
                      help="Output both bases bedGraph (stdout) and windows"
                           "bed12 (specified file).")
    parser.add_option("-t", "--threshold", dest="threshold", type="float",
                      default=0.05,
                      help="p-value threshold under which to merge windows")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # Standard in contains the transcripts
    
    gffs = GTF.gene_iterator(GTF.iterator(options.stdin))

    # bam file is the first positional arguement
    bamfile = pysam.Samfile(args[0])

    if options.output_both:
        outfile_bases = options.stdout
        outfile_windows = IOTools.openFile(options.output_both, "w")
    elif options.output_windows:
        outfile_bases = None
        outfile_windows = options.stdout
    else:
        outfile_bases = options.stdout
        outfile_windows = None

    if options.fdr and options.pipeout:
        E.warning("--fdr implies not --pipeout, instant output disabled")
        options.pipeout = False

    if options.pipeout:
        output = InstantOutput(outfile_bases=outfile_bases,
                               outfile_windows=outfile_windows,
                               window_size=options.window_size,
                               threshold=options.threshold)
    else:
        output = DeferredOutput(outfile_bases=outfile_bases,
                                outfile_windows=outfile_windows,
                                correct=options.fdr,
                                window_size=options.window_size,
                                threshold=options.threshold)

    E.info("Counting accross transcripts ...")
    max_end = 0
    for gene in gffs:

        if options.grouping == "all":
            gene = GTF.merged_gene_iterator(gene)

        transcript_ps = {}

        for transcript in gene:
            
            # E.debug("Transcript is %s" % transcript[0].transcript_id)
            coords_converter = iCLIP.TranscriptCoordInterconverter(transcript)
            exons = GTF.asRanges(transcript, "exon")
            counts = iCLIP.count_intervals(bamfile,
                                           exons,
                                           strand=transcript[0].strand,
                                           contig=transcript[0].contig,
                                           dtype=options.dtype)

 
            counts.index = coords_converter.genome2transcript(counts.index.values)
            counts = counts.sort_index()
            cds = GTF.asRanges(transcript, "CDS")

            if options.grouping == "utrs" and len(cds) > 0:
                
                cds_interval = (cds[0][0], cds[-1][1])
                cds_interval = coords_converter.genome2transcript(cds_interval)
                cds_interval.sort()
                cds_length = cds_interval[1] - cds_interval[0]

                p_intervals = [(0, cds_interval[0]),
                               (cds_interval[0], cds_length),
                               (cds_interval[1], coords_converter.length - cds_interval[1])]

            else:  # do not group by cds or there is no cds
                p_intervals = [(0, coords_converter.length)]

            p_values = [calculateProbabilities(counts, options.window_size,
                                              length=length, start=start)
                        for start, length in p_intervals
                        if length > 0]
  
            if len(p_values) > 1:
                p_values = pd.concat(p_values)
            else:
                p_values = p_values[0]

            p_values.index = coords_converter.transcript2genome(p_values.index.values)
 
 
            intron_intervals = GTF.toIntronIntervals(transcript)
            
            if len(intron_intervals) > 0:
                intron_coords = iCLIP.TranscriptCoordInterconverter(transcript,
                                                                    introns=True)
                intron_counts = iCLIP.count_intervals(bamfile,
                                                      intron_intervals,
                                                      strand=transcript[0].strand,
                                                      contig=transcript[0].contig,
                                                      dtype=options.dtype)
             
                intron_counts.index = intron_coords.genome2transcript(
                    intron_counts.index.values)
                intron_counts = intron_counts.sort_index()
                intron_pvalues = calculateProbabilities(intron_counts,
                                                        options.window_size,
                                                        intron_coords.length)
                                                        
                intron_pvalues.index = intron_coords.transcript2genome(
                    intron_pvalues.index.values)
                p_values = p_values.append(intron_pvalues)
                
            transcript_ps[transcript[0].transcript_id] = p_values

        transcript_df = pd.DataFrame(transcript_ps)

        transcript_df.index.rename("position", inplace=True)
        transcript_df["contig"] = gene[0][0].contig
        transcript_df["strand"] = gene[0][0].strand
        transcript_df["gene_id"] = gene[0][0].gene_id
        transcript_df.set_index("contig", append=True, inplace=True)
        transcript_df.set_index("strand", append=True, inplace=True)
        transcript_df.set_index("gene_id", append=True, inplace=True)
 
        gene_ps = transcript_df.mean(1)
        gene_ps = gene_ps.reorder_levels(["gene_id", "contig",
                                          "strand", "position"])

        output.write(gene_ps, gene)

    output.close()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
