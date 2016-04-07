'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script calculates the "reproducibility" for a set of iCLIP BAM files.
The reproducibility at level n for file i is the fraction of bases with 
coverage n in file i that have at least coverage 1 in one of the other files.

Memory and CPU usage should be reasonably small and grow in proportion to a) The 
number of samples and b) The number of unique bases with reads mapped (
and not the number of reads in the input file).

Options
-------

--dtype, Count data is stored internally in this data type. The default is
        uint16 That is 16bit unsigned integers. Using a smaller integer type,
        like uint8 might save memory, but limits the largest possible count
        on any one base.

-m, --max-level, This allows the specification of what "level" of reproducibilty 
        to calculate up to. By default, this is the maximum depth found for a
        particular track. However, the counts at very high depths are probably
        not of a great amount of use, and so the max depth can be limited to
        save time.

-c, --contig, This allows the restricting of the calculation to a single 
        chromosome. This could be useful if it was neccesary to parrellise the
        excution for any reason, or for quick testing perpuses. 

Usage
-----
<Example use case>

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pysam
import iCLIP
import CGAT.Experiment as E
import numpy as np
import collections
import pandas as pd
import os.path


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--dtype", dest="dtype", type="string",
                      help="Numpy dtype for storing counts in. [%default]",
                      default="uint16")
    parser.add_option("-m", "--max-level", dest="max_level",
                      help="Max level to calculate reproducibility at."
                           " Use 0 for max possible for each track (default)",
                      default=0)
    parser.add_option("-c", "--contig", dest="contig",
                      help="Restrict analysis to a particular chromosome",
                      default=None)
    parser.add_option("-t", "--track", dest="track",
                       help="Restrict analysis to one of the input samples vs."
                            "all the others",
                       default=None)
        
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    samfiles = [pysam.Samfile(fn, 'rb') for fn in args]
#    total_counter = [E.Counter() for samfile in samfiles]
    running_totals = {sf: [collections.defaultdict(int)
                           for x in range(len(args) - 1)]
                      for sf in args}
    running_hits = {sf: [collections.defaultdict(int)
                         for x in range(len(args) - 1)]
                    for sf in args}

    contigs = zip(samfiles[0].references, samfiles[0].lengths)
    
    if options.track:
        use_index, use_names = zip(*[(i,fn) for i,fn in enumerate(args)
                               if fn == options.track])
    else:
        use_index, use_names = zip(*enumerate(args))

    if options.contig:
        contigs = [x for x in contigs if x[0] == options.contig]

    E.debug("Reporting on input(s) %s: %s" %(",".join(map(str,use_index)),",".join(use_names)))

    for ref, length in contigs:

        E.debug("Starting %s, length %i" % (ref, length))

        depths = pd.DataFrame(dtype=options.dtype)

        for sf in range(len(samfiles)):
            E.debug("Reading File %s" % args[sf])
            pos_depth, neg_depth, counter = \
                iCLIP.countChr(samfiles[sf].fetch(ref), length, options.dtype)
            #(depths[:length, sf], depths[length:, sf], counter) = \
            #    iCLIP.countChr(samfiles[sf].fetch(ref), length, options.dtype)
            
            neg_depth.index = neg_depth.index + length
            depth = pd.concat([pos_depth, neg_depth])
            depth.name = args[sf]
            
            depths = depths.join(depth, how="outer")
#            total_counter[sf] += counter
        
        depths = depths.fillna(0)

        for sf in use_index:
        
            try:
                E.debug("Max depth for %s is %i" % 
                        (args[sf], int(depths.iloc[:, sf].max())))
            except ValueError:
                E.warn("Zero max depth for both")
                continue



            if int(options.max_level) == 0:
                n_max = int(depths.iloc[:, sf].max())
            else:
                n_max = int(options.max_level)

            for n in range(n_max):

                E.debug("Calculating %i level reproducibility for file %s"
                        % (n, args[sf]))

                sites = depths.iloc[:, sf] > n
                for i in range(len(args) - 1):
                    running_totals[args[sf]][i][n] += sites.sum()
                
                replicating_sites = \
                    depths.ix[sites, np.arange(len(samfiles)) != sf] > 0
                n_replicating_samples = replicating_sites.sum(axis=1)

                for i in range(len(args) - 1):
                    running_hits[args[sf]][i][n] += \
                        (n_replicating_samples > i).sum()
                del sites
                del replicating_sites
                del n_replicating_samples

        del depths

    outlines = []
    for sf in use_names:
        for i in range(len(running_totals[sf])):
            for n in running_totals[sf][i]:
                outlines.append([os.path.basename(sf), i+1, n+1, running_hits[sf][i][n],
                                 running_totals[sf][i][n]])

    header = ["Track", "fold", "level", "hits", "totals"]

    outlines = "\n".join(["\t".join(map(str, line)) for line in outlines])
    outlines = "\t".join(header) + "\n" + outlines + "\n"

    options.stdout.write(outlines)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
