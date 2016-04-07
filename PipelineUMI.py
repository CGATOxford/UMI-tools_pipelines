import iCLIP
import pandas
import pysam

import CGAT.Experiment as E
from CGAT import IOTools
from CGATPipelines.Pipeline import cluster_runnable
from CGAT import Bed

@cluster_runnable
def getSigHeights(sig_bed, bam_file, outfile):
    ''' Take a bedgraph of significant x-linked bases and return a begraph of heights '''

    print iCLIP.__file__
    bam = pysam.AlignmentFile(bam_file)
    last_contig = None
    intervals = []
    outlist = []
    for line in IOTools.openFile(sig_bed):

        contig, start, end, pval = line.strip().split("\t")

        if last_contig is None:
            last_contig = contig

        if contig != last_contig:
            print last_contig
            out = iCLIP.count_intervals(bam, intervals, last_contig)
            out.index.name = "start"
            out.name = "count"
            out = out.reset_index()
            out["contig"] = last_contig
            outlist.append(out)
            intervals = []
            last_contig = contig
        
        intervals.append((int(start), int(end)))


    # output the final chrom
    out = iCLIP.count_intervals(bam, intervals, last_contig)
    out.index.name = "start"
    out.name = "count"
    out = out.reset_index()
    out["contig"] = last_contig
    outlist.append(out)
    
    outframe = pandas.concat(outlist)
    outframe["end"] = outframe["start"] + 1
    outframe = outframe.loc[:,["contig", "start", "end", "count"]]
    outframe["start"] = outframe["start"].astype("int")
    outframe["end"] = outframe["end"].astype("int")

    outframe.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False, header=False)

    
@cluster_runnable
def countTagsInClusters(bedfile, bamfile, outfile):

    bam = pysam.AlignmentFile(bamfile)

    outlines = []

    for bed in Bed.iterator(IOTools.openFile(bedfile)):
        interval = (bed.start, bed.end)
        counts = iCLIP.count_intervals(bam, [interval], bed.contig).sum()
        outlines.append(["%s:%i-%i" % (bed.contig, bed.start, bed.end), str(counts)])

    IOTools.writeLines(outfile, outlines, header=["position","count"])
