import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGAT.FastaIterator as FastaIterator
import CGAT.Experiment as E
#import CGATPipelines.PipelineUtilities as PUtils
from CGATPipelines.Pipeline import cluster_runnable
import pandas
import os
import re
import pysam

# The PARAMS dictionary must be provided by the importing
# code

PARAMS = {}


def checkParams():

    if not len(PARAMS) > 0:
        raise ValueError(
                "Please set PARAMS dictionary in PipelineiCLIP module")


def removeFirstAndLastExon(infile, outfile):

    transcripts = GTF.transcript_iterator(
        GTF.iterator(IOTools.openFile(infile)))
    outfile = IOTools.openFile(outfile, "w")

    for transcript in transcripts:

        for exon in transcript[1:-1]:
            outfile.write(str(exon) + "\n")

    outfile.close()


def getBarcodeCG(table, outfile):
    ''' Annotate barcode use statistics with %GC '''

    statement = " SELECT * FROM %(table)s" % locals()

    umi_stats = PUtils.fetch_DataFrame(statement)

    def _GC(x):
        return float(x.count("G") + x.count("G"))/len(x)

    barcode_gc = umi_stats.Barcode.apply(_GC)
    sample_gc = umi_stats.Sample.apply(_GC)
    umi_gc = umi_stats.UMI.apply(_GC)

    gc_stats= pandas.DataFrame({"Barcode":umi_stats.Barcode,
                                "barcode_gc": barcode_gc,
                                "sample_gc": sample_gc,
                                "umi_gc": umi_gc})

    gc_stats.to_csv(IOTools.openFile(outfile,"w"), sep="\t", index=False)


###################################################################
def callClusters(bamfile, gtffile, outfiles,
                 window_size=None,
                 pthresh=None):
    ''' Wrapper around find_reproducible_clusters.
    If no window_size is specified, it is taken from
    pipeline.ini'''

    checkParams()

    bedGraph, bed12 = outfiles
    logfile = P.snip(bed12, ".bed.gz")

    if window_size:
        options = "--window-size=%i" % window_size
    else:
        options = "--window-size=%s" % PARAMS["clusters_window_size"]

    if PARAMS["clusters_fdr"]:
        options += " --fdr"
    if PARAMS["clusters_grouping"]:
        options += " --grouping=%s" % PARAMS["clusters_grouping"]

    if pthresh:
        options += " -t %s" % str(pthresh)
    else:
        options += " -t %s" % PARAMS["clusters_pthresh"]

    
    job_options = "-l mem_free=10G"
    statement = '''python %(scriptsdir)s/gtf2gtf.py -L %(logfile)s.log
                           -I %(gtffile)s
                          --method=sort --sort-order=gene+transcript
                 | python %(scriptsdir)s/gtf2gtf.py -L %(logfile)s.log 
                          --method=set-transcript-to-gene
                 | python %(project_src)s/find_significant_bases.py
                   %(bamfile)s
                   %(options)s
                   --output-both=%(bed12)s
                  -L %(logfile)s.log
                | gzip -c > %(bedGraph)s '''

    P.run()


###################################################################
def callReproducibleClusters(infiles, outfile, min_overlap):
    '''Find clusters that appear in more than one replicate'''

    checkParams()

    merge_template = '''<( zcat %s
                          | sort -k1,1 -k2,2n
                          | python %s/bed2bed.py
                          --method=merge
                          --merge-and-resolve-blocks
                          --merge-stranded
                           -L /dev/null
                            ) '''
    infiles = " ".join(
        [merge_template % (infile, PARAMS["scriptsdir"])
         for infile in infiles])
    logfile = P.snip(outfile, ".bed.gz")
    statement = ''' cat %(infiles)s
                  | sort -k1,1 -k2,2n -k3,3n
                  | python %(scriptsdir)s/bed2bed.py
                          --method=merge
                          --merge-and-resolve-blocks
                          --merge-min-intervals=%(min_overlap)s
                          --merge-stranded
                           -L %(logfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
def removeInputOverlappingClusters(sample, control, outfile):
    '''Remove reproducible clusters that overlap with reproducible
    input clusters '''

    statement = ''' bedtools intersect -v -a %(sample)s -b %(control)s
                   | gzip > %(outfile)s '''
    P.run()


###################################################################
def clustersToBigBed(infile, outfile):
    '''Convert beds to bigbed '''

    checkParams()

    tmp = P.getTempFilename()
    genome_file = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_contigs_tsv"])
    statement = ''' zcat %(infile)s | sort -k1,1 -k2,2n 
                    | awk 'BEGIN{OFS="\\t"} $5=1' > %(tmp)s;
                    checkpoint;
                    bedToBigBed %(tmp)s %(genome_file)s %(outfile)s;
                    checkpoint;
                    rm %(tmp)s'''
    P.run()


###################################################################
def makeClustersUCSC(infiles, outfile, group, label):
    '''Compile UCSC track file for bigbeds from list of cluster files'''

    template = '''
         track %(track_name)s
         parent %(group)s %(visible)s
         shortLabel %(short_label)s
         longLabel %(long_label)s
         bigDataUrl %(big_data_url)s
         type bigBed 12'''

    outlines = []

    for infile in infiles:

        big_data_url = os.path.basename(infile)
        

        if "reproducible" in infile:
            visible = "on"
            track_name = group + "_" + re.match(
                ".+/(.+).reproducible.*bigBed", infile).groups()[0]
            long_label = "Clusters from %s appearing in at least %s replicates" \
                         % (track_name, PARAMS["clusters_min_reproducible"])
        else:
            visible = "off"
            track_name = group + "_" + re.match(
                ".+/(.+)\.bigBed", infile).groups()[0]
            long_label = "Clusters from track %s" % track_name

        short_label = "%s clusters" % track_name
        outlines.append(template % locals())

    composite_stanaz = '''

    track %(group)s
    superTrack on
    shortLabel %(label)s
    longLabel %(label)s
    '''
    outlines = [composite_stanaz % locals()] + outlines

    outlines = "\n".join(outlines)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(outlines+"\n\n")


###################################################################
def subsampleNReadsFromFasta(infile, outfile, nreads, logfile=""):

    checkParams()

    nseqs = FastaIterator.count(infile)

    if nreads > nseqs:
        prop = 1
    else:    
        prop = float(nreads)/float(nseqs)

    if logfile:
        logfile = "-L %s" % logfile

    statement = ''' python %(scriptsdir)s/fasta2fasta.py 
                     -I %(infile)s
                     %(logfile)s
                     -m sample
                     --sample-proportion=%(prop)s
                     -S %(outfile)s '''

    P.run()


###################################################################
@cluster_runnable
def calculateSplicingIndex(bamfile, gtffile, outfile):

    bamfile = pysam.AlignmentFile(bamfile)

    counts = E.Counter()

    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(gtffile))):

        introns = GTF.toIntronIntervals(transcript)
        E.debug("Gene %s (%s), Transcript: %s, %i introns" %
                (transcript[0].gene_id,
                 transcript[0].contig,
                 transcript[0].transcript_id,
                 len(introns)))

        for intron in introns:
            reads = bamfile.fetch(
                reference=transcript[0].contig,
                start=intron[0], end=intron[1])
            
            for read in reads:
                if 'N' in read.cigarstring:
                    blocks = read.get_blocks()
                    starts, ends = zip(*blocks)
                    if intron[0] in ends and intron[1] in starts:
                        counts["Exon_Exon"] += 1
                    else:
                        counts["spliced_uncounted"] += 1
                elif (read.reference_start <= intron[0] - 3
                      and read.reference_end >= intron[0] + 3):
                    if transcript[0].strand == "+":
                        counts["Exon_Intron"] += 1
                    else:
                        counts["Intron_Exon"] += 1
                elif (read.reference_start <= intron[1] - 3
                      and read.reference_end >= intron[1] + 3):
                    if transcript[0].strand == "+":
                        counts["Intron_Exon"] += 1
                    else:
                        counts["Exon_Intron"] += 1
                else:
                    counts["unspliced_uncounted"] += 1

        E.debug("Done, counts are: " + str(counts))
    header = ["Exon_Exon",
              "Exon_Intron",
              "Intron_Exon",
              "spliced_uncounted",
              "unspliced_uncounted"]

    with IOTools.openFile(outfile, "w") as outf:

        outf.write("\t".join(header)+"\n")
        outf.write("\t".join(map(str, [counts[col] for col in header]))
                   + "\n")


        

