'''
iCLIP_bam2geneprofile.py - produce geneprofile of iCLIP sites
============================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

While bam2geneprofile in CGAT is a very flexible tool, it is not
neccesarily suitable for iCLIP as we should only consider the first
base (or any mutant bases).

This script wraps iCLIP.meta_gene to produce metagene profiles of
iCLIP bam files. In future it may offer ability to work from bigwigs
or beds.

Using single bases means we don't have to worry about over sampling single
reads, so profiles should be less resolution sensitive

Usage
-----

.. Example use case

Example::

   python iCLIP_bam2geneprofile.py -I geneset.gtf.gz mybam.bam

Type::

   python iCLIP_bam2geneprofile.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import iCLIP
import pysam
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--output-matrix", dest="matrix", type="string",
                      default=None,
                      help="output full matrix to this file")
    parser.add_option("-f", "--flanks", dest="flanks", type="int",
                      default=100,
                      help="number of basepairs to use for gene flanks")
    parser.add_option("-b", "--exon-bins", dest="exon_bins", type="int",
                      default=1000,
                      help="number of bins to divide transcripts into")
    parser.add_option("--flank-bins", dest="flank_bins", type="int",
                      default=10,
                      help="number of bins to divide flanks into")
    parser.add_option("--scale-flanks", dest="scale_flanks", action="store_true",
                      default=False,
                      help="Scale the size of the flank bins to match the size of the"
                      "exon bins for each transcript")
    parser.add_option("--pseudo_count", dest="pseudo_count", type="float",
                      default=0,
                      help="add pseduo count to bins to mitiage effects of low numbers of reads")
    parser.add_option("--normalised_profile", dest="normalize_profile", action="store_true",
                      default=False,
                      help="Normlize profile by profile sum")
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    bam = pysam.AlignmentFile(args[0])
    
    if options.flanks > 0:
        bins = [options.flank_bins,
                options.exon_bins,
                options.flank_bins]
    else:
        bins = options.exon_bins

    summed_matrix, counts_matrix = iCLIP.meta_gene(
        options.stdin,
        bam,
        bins,
        options.flanks, 
        output_matrix=options.matrix is not None,
        calculate_flanks=options.scale_flanks,
        pseudo_count=options.pseudo_count)

    print summed_matrix
    try:
        summed_matrix = summed_matrix[["flank5", "exons", "flank3"]]
    except IndexError:
        pass

    summed_matrix = summed_matrix.reset_index()

    if options.normalize_profile:
        summed_matrix["density"] = summed_matrix["density"]/summed_matrix["density"].sum()

    summed_matrix.to_csv(options.stdout, sep="\t",
                         index=True,
                         index_label="bin")

    if options.matrix:
        counts_matrix = counts_matrix.transpose()
        print counts_matrix.index
        counts_matrix = counts_matrix.loc[["flank5", "exons", "flank3"],:]
        counts_matrix = counts_matrix.reset_index(drop=True)
        counts_matrix = counts_matrix.transpose()

        counts_matrix.to_csv(IOTools.openFile(options.matrix, "w"),
                             sep="\t",
                             index=True,
                             index_label="transcript_id")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
