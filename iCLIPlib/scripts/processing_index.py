'''
processing_index.py - calculate 3' processing index on an iCLIP BAM
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Wraps meta.processing_index. Calculates the "Processing index" for the supplied
BAM files. The procesing index is a measure of whether transcripts are bound
to the primary or mature mRNA, and is defined as:

.. math::
   pi = log_2(\frac{\sum_{i=1}^G N^{PM}}{\sum_{i=1}^G N^M})

where :math:`N^{PM}` is the number tags on the pre-mRNA and :math:`N^M` is the
tags on mature mRNA and G is the number of genes in the set. 

These are calculated by counting tags in a window (default 50bp) upstream and downstream 
of a cleavage site. Because upstream of the cleavage site is a mix of mature
and pre-mRNA, while downstream is only pre-mRNA, then :math:`N^{PM} == N^{down} - N^{up}`

This is after Baejen et al Mol Cell 5(55):745-757, however not quite the same
as they include a normalisation factor of 1/G which I don't think makes sense.


Usage
-----

Bed or GTF formatted cleavage site annotations come in from stdin. BAM file is
specified as the positional arguement

.. Example use case

Example::

   python processing_index.py -I BED/GTF/GFF [OPTIONS] BAMFILE

Type::

   python processing_index.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pysam
import os

import CGAT.Experiment as E
from CGAT import GTF
from CGAT import Bed

sys.path.insert(1, os.path.join(
    os.path.dirname(__file__), ".."))

import iCLIP


def last_exon_transcript(gff_file):
    for transcript in GTF.transcript_iterator(GTF.iterator(gff_file)):
        transcript = sorted(transcript)
        if transcript[0].strand == "-":
            yield transcript[0]
        else:
            yield transcript[-1]


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-f", "--input-format", type="choice",
                      choices=["bed", "gtf", "gff"],
                      default="bed",
                      dest="format",
                      help="Format of the cleavage site definition"
                           " default [%default]")
    parser.add_option("--feature", type="choice",
                      choices=["gene", "transcript", "entry"],
                      default="gene",
                      dest="feature",
                      help="Which feature to use if using gtf")
    parser.add_option("-w", "--window-size", type="int",
                      default=50,
                      dest="window_size",
                      help="Number of bases to count upstream and downstream"
                           " of cleavage site. [%default]")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    try:
        bamfile = pysam.AlignmentFile(args[0])
    except IndexError:
        E.error("Please supply a bam file as the first positional arguement")
        return 1
    except IOError:
        E.error("Cannot open BAM file %s" % args[0])

    interval_iterators = {"bed": Bed.iterator,
                          "gtf-gene": lambda x:
                          GTF.merged_gene_iterator(GTF.iterator(x)),
                          "gtf-transcript": last_exon_transcript}

    if options.format == "gtf":
        options.format += "-" + options.feature

    iterator = interval_iterators[options.format](options.stdin)

    pi = iCLIP.processing_index(iterator, bamfile, options.window_size)

    options.stdout.write("Processing Index\t%s\t%s\n" % (args[0], pi))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
