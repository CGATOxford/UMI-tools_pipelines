'''
reproducbility_by_exon.py - calculate reproducbility per exon
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Use various methods to calculate the similarity between two tracks. In each
case the metric is calculated for the exon. The significance of this is 
calculated by performing randomisations of the first profile. The standard
deviations of the randomisaitons will be calculated and the Z score returned.

Average Distance
++++++++++++++++

The distance between all pairwise combinations of reads are found and the average
taken. 

Min Distance
++++++++++++

For each tag in profile one, the closest tag in profile two is found and the
average taken. 

Correlation
+++++++++++

The tags in each profile are expanded (default 10 bases) and the spearman 
correlation between the resulting profiles calculated. 

Usage
-----

.. Example use case

Example::

   python reproducbility_by_exon.py PROFILE_1_FILE PROFILE_2_FILE -I GTF_FILE

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import pysam
import iCLIP

def bed_counter(bedfile, intervals, contig, strand):

    intervals = []

    for interval in intervals:
        tmp_intervals = bedfiles[contig].fetch(*interval)
	intervals.extend([tint[0] for tint in tmp_intervals
                        if tint[0] >= interval[0] and
                           tint[1] < interval[1] and
                           tint[2].strand == strand])

    profile = pd.Series(intervals)
    profile = profile.count_values()
    profile = profile.sort_index()
    return profile

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv


    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m","--method", dest="method", type="choice",
                      choices=["ave_dist","min_dist","corr"],
                      default="min_dist",
                      help="Method for calcuating similarity between profiles")
    parser.add_option("-s", "--spread", dest="spread", type="int",
                      default=10,
                      help="Amount to spread each tag by")
    parser.add_option("-k", "--keep-dist", dest="keep_dist", 
                      action="store_true",
                      help="Keep the distribution of tag depths")
    parser.add_option("-r", "--rands", dest="rands", 
                      default=100,
                      help="Number of randomisations to use for calculating"
                           " mean and stdev of distance")
 
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    profile1_file, profile2_file = args
    profile1_file = pysam.AlignmentFile(profile1_file)

    if profile2_file.endswith("bed") or profile2_file.endswith("bed.gz"):
	profile2_file = Bed.readAndIndex(profile2_file, with_values=True)
        profile2_counter = bed_counter
    else:
        profile2_file = pysam.AlignmentFile(profile2_file)
        profile2_counter = iCLIP.count_intervals

    if options.method=="min_dist":
        distance_func = iCLIP.findMinDistance
    elif options.method=="ave_dist":
        distance_func = iCLIP.calcAverageDistance
    else:
        def distance_func(profile1,profile2):
            return iCLIP.corr_profile(profile1,profile2, options.spread, profile2_ready=True)
 
    for exon in GTF.iterator(options.stdin):
	if exon.feature != "exon":
            continue

        contig = exon.contig
        strand = exon.strand
        transcript_id = exon.transcript_id
        start = exon.start
        end = exon.end

        profile1 = iCLIP.count_intervals(profile1_file, 
                                         [(start, end)],
                                         contig=contig,
                                         strand=strand)

        profile2 = profile2_counter(profile2_file,
                                    [(start, end)],
                                    contig=contig,
                                    strand=strand)

        if profile1.sum() == 0 or profile2.sum() == 0:
            z = "NA"
            distance = "NA"
            options.stdout.write(
                "%(contig)s\t%(start)i\t%(end)i\t%(transcript_id)s\t%(strand)s\t%(distance)s\t%(z)s\n" % locals())
            continue

	if options.method=="corr":
             profile2 = iCLIP.spread(profile2, options.spread)

        distance = distance_func(profile1, profile2)

        rands = iCLIP.rand_apply(profile=profile1, 
                                 exon=exon, 
                                 n=options.rands, 
                                 func=distance_func,
                                 keep_dist=options.keep_dist,
                                 profile2=profile2)

        z = (distance - rands.mean())/rands.std()

        options.stdout.write(
          "%(contig)s\t%(start)i\t%(end)i\t%(transcript_id)s\t%(strand)s\t%(distance).3f\t%(z).2f\n" % locals())
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
