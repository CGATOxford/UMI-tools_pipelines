'''
iCLIP2bigWig -- convert iCLIP BAM files to two wig files
============================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Takes a bam file on the stdin and converts it into two
bigWig files, each containing the depth of crosslinked bases
at any given position on each strand.

Output files are named according to the provided template
with _plus and _minus suffixes.

Requires the ucsc wigToBigWig program to work properly.

Options
-------

Will read bamfile off of the stdin or from -I, but 
BAM file must be indexed and so cannot be manipulated on 
stdin.

If wig files are required as output, for example if 
wigToBigWig is not installed, --wig will output wig files
rather than bigWig files.


Usage
-----


python umi_stats.py -I [BAMFILE] [OUT_TEMPLATE]



Command line options
--------------------

'''

import sys
import os
import shutil
import CGAT.Experiment as E
import pysam
import subprocess
import tempfile
import iCLIP


def outputToBW(infile, outfile_prefix, chrom_sizes):

    E.debug("Attempting to output %s to bigWig" % outfile_prefix)
    command = ["wigToBigWig",
               infile,
               chrom_sizes,
               outfile_prefix + ".bw"]
    try:
        subprocess.check_call(command)
    except Exception as e:
        E.error("Error on conversion: %s" % e)
        E.info("Outputting to wig")
        shutil.move(infile, outfile_prefix + ".wig")
    else:
        E.debug("Conversion successful")
        os.unlink(infile)
    
def outputToWig(depths,chrom, wigfile):
    '''depths is a pandas series keyed on chromosome position,
    chrom is a chromosome, wigfile is a file to output to.
    This function converts a series of depths into wig formated
    text and writes it to the specified file '''

    wigfile.write("variableStep\tchrom=%s\n" % chrom)
    for row in depths.iteritems():
        row = list(row)
        row = "\t".join(map(str,row)) + "\n"
        wigfile.write(row)

def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-w","--wig", dest="output_wig", action="store_true",
                      default=False,
                      help="Write output to wig file rather than bigwig")
    parser.add_option("--dtype", dest = "dtype", type="string",
                      default="uint32",
                      help="dtype for storing depths")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.stdin == sys.stdin:
        in_bam = pysam.Samfile("-", "rb")
    else:
        fn = options.stdin.name
        options.stdin.close()
        in_bam = pysam.Samfile(fn, "rb")

    plus_wig = tempfile.NamedTemporaryFile(delete=False)
    minus_wig = tempfile.NamedTemporaryFile(delete=False)

    contig_sizes = []

    for chrom, chrom_length in zip(in_bam.references, in_bam.lengths):

        # get depths over chromosome
        pos_depth, neg_depth, counter = iCLIP.countChr(in_bam.fetch(chrom),
                                                       chrom_length,
                                                       options.dtype)
        pos_depth_sorted = pos_depth.sort_index()
        del pos_depth
        neg_depth_sorted = neg_depth.sort_index()
        del neg_depth
        neg_depth_sorted = -1*neg_depth_sorted

#        E.debug("Counted %i truncated on positive strand, %i on negative"
#                % (counter.truncated_pos, counter.truncated_neg))
#        E.debug("and %i deletion reads on positive strand, %i on negative"
#                % (counter.deletion_pos, counter.deletion_neg))

        # output to temporary wig file
        outputToWig(pos_depth_sorted, chrom, plus_wig)
        outputToWig(neg_depth_sorted, chrom, minus_wig)
    
        contig_sizes.append([chrom, chrom_length])

        del pos_depth_sorted
        del neg_depth_sorted
    
    plus_wig_name = plus_wig.name
    minus_wig_name = minus_wig.name
    plus_wig.close()
    minus_wig.close()

    outname_plus = args[0] + "_plus"
    outname_minus = args[0] + "_minus"

    if options.output_wig:
        E.debug("Outputting to wig")
        shutil.move(plus_wig_name, outname_plus + ".wig")
        shutil.move(minus_wig_name, outname_minus + ".wig")
    else:
        chrom_sizes_file = tempfile.NamedTemporaryFile(delete=False)
        contig_sizes = ["\t".join(map(str,row)) for row in contig_sizes]
        contig_sizes = "\n".join(contig_sizes) + "\n"
        chrom_sizes_file.write(contig_sizes)
        chrom_sizes_filename = chrom_sizes_file.name
        chrom_sizes_file.close()

        outputToBW(plus_wig_name, outname_plus, chrom_sizes_filename)
        outputToBW(minus_wig_name, outname_minus, chrom_sizes_filename)


    # write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
