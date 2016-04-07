'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import CGAT.Experiment as E
import CGAT.Bed as Bed
import pandas
import iCLIP
import pysam
import numpy


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")
    parser.add_option("-w", "--window-size", dest="window_size", type="int",
                      help="size of window to summarise over")
    parser.add_option("-c", "--contig", dest="contig", type="string",
                      help="restrict to contig")
    parser.add_option("--dtype", dest="dtype", default = "uint32")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    try:
        samfile = pysam.Samfile(args[0], "rb")
    except IndexError:
        raise ValueError("Please supply a BAM file as the first arguement")

    contigs = zip(samfile.references, samfile.lengths)
    if options.contig:
        contigs = [x for x in contigs if x[0] == options.contig]

    for contig, length in contigs:

        E.debug("Doing chromosome %s ..." % contig)
        E.debug("Getting depth vector")
        pos_depths, neg_depths, _ = \
            iCLIP.countChr(samfile.fetch(contig), length, options.dtype)
        
        E.debug("Binning counts ...")
        num_bins = -(-length//options.window_size)
        bin_edges = numpy.arange(num_bins+1) * options.window_size
        
        pos_bin_sums = pos_depths.groupby(
            pos_depths.index.values//options.window_size).sum()
        neg_bin_sums = neg_depths.groupby(
            neg_depths.index.values//options.window_size).sum()

        E.debug("Creating bed entries ...")
        def _score2bed(bin, score, strand):

            start = int(bin)*options.window_size
            end = start + options.window_size
            return "\t".join([contig, str(start), str(end),
                              ".", str(int(score)), strand])


        for bin, score in pos_bin_sums.iteritems():
            options.stdout.write(_score2bed(bin, score, "+") + "\n")

        for bin, score in neg_bin_sums.iteritems():
            options.stdout.write(_score2bed(bin, score, "-") + "\n")
        
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
