'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Test the enirhcment on kmers in regions surrounding iCLIP sites, using
the methodology described in Wang et al, 2010, PLoS Biol. 8:e1000530 (
:PMID:`21048981`).


Usage
-----

   python iCLIP_kmer_enrichment.py -b BAMFILE -f FASTA -I GTFFILE [OPTIONS]

Example::

   zcat gene_models.gtf.gz | python iCLIP_kmer_enrichment.py -b my_bam.bam -f genome.fa.gz

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
from CGAT.IndexedFasta import IndexedFasta
import pysam

import iCLIP

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--bam-file", dest="bam", type="string",
                      help="BAM file with iCLIP reads")
    parser.add_option("-f", "--fasta-file", dest="fasta", type="string",
                      help="CGAT indexed Fasta file with genome sequence")
    parser.add_option("-k", "--kmer", dest="kmer", type="int",
                      default=5,
                      help="Size of kmer to test default=[%default]")
    parser.add_option("-s", "--spread", dest="spread", type="int",
                      default=15,
                      help="Amount of sequence around each read to consider"
                      "default=[%default]")
    parser.add_option("-n", "--num-randomizations", dest="randomisations",
                       type="int", default=100,
                       help="Number of times to permute profiles to assess"
                            "significance of enrichment")
    parser.add_option("-p", "--processes", dest="proc", type="int",
                      default=None,
                      help="Use this many processesors for multiprocessing")
    parser.add_option("--feature", dest="feature", type="choice",
                      choices=["transcript", "gene"],
                      default="gene",
                      help="Treat transcripts seperately or merge genes")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.proc:
        try:
            import multiprocessing as mp
            pool = mp.Pool(options.proc)
        except ImportError:
            E.warn("Multiprocessing setup failed."
                   " Falling back to single processor mode")
            pool = None
    else:
        pool = None

    bam = pysam.AlignmentFile(options.bam)
    fasta = IndexedFasta(options.fasta)

    if options.feature == "gene":
        gtf_iterator = GTF.flat_gene_iterator(
            GTF.iterator(options.stdin))
    else:
        gtf_iterator = GTF.transcript_iterator(
            GTF.iterator(options.stdin))

    results = iCLIP.pentamer_enrichment(gtf_iterator,
                                        bam,
                                        fasta,
                                        options.kmer,
                                        options.randomisations,
                                        spread=options.spread,
                                        pool=pool)

    results.name = "Z"

    results.to_csv(options.stdout, header=True, index_label="Kmer", sep="\t")
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
