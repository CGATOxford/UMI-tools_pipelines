###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""
===========================
Pipeline iCLIP
===========================

:Author: Ian Sudbery
:Release: 0.0.1
:Date: 18/12/2016
:Tags: Python



Overview
========
This pipeline executes analysis to produce the results in 
in Smith *et al* 2017 that deal with iCLIP data. The pipeline runs in a simiar 
way to a make pipeline. Once run, the pipeline will
automatically create a succession of jobs that run the required analysis tasks.
At each step the pipeline tests if the input files for a task exist, if do, the
date stamps on the output files are compared to the input files and if the input
flies are newer, the task is run to regenerate the output files. If the input
files are not present, they task that generates them is activated. 

Usage
=====

Initialize the pipeline by creating a directory with the input files (see input below) 
and running::
    python PATH_TO_SRC/pipeline_iCLIP.py config

This will copy the required configuration files into the directory (see configuration)
below. Once setup and configured (see below) the pipeline is executed
with::
    python PATH_TO_SRC/pipeline_iCLIP.py make full

By default the pipeline will attempt to use drmaa to execute pipeline tasks
on a cluster. The default cluster manager is SGE, but SLRUM is also possible. 
Cluster execution is the most efficient way to run the pipeline, and on our
system, running 100 cluster jobs, allocating 12 cores to each bowtie run and 6 
cores to each peakcalling run the pipeline takes 50 hours. It is also 
possible to run without a cluster by specifying the `--no-cluster` option. The 
number of concurent jobs to be run can be set with `-p`. Running a single 
job on a local machine, we estimate that the pipeline could take several 
months to run, and would require at least 32GB of RAM. 

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The default
values will produce the analyses in the paper for the most part, but some
site specific values must be set:

* The location of the bowtie or star index files and reference genome
* The location of a GTF file with the relevant gene annotation. The files used
  for the SRSF analysis are included in the data directory of the repository,
  edit the configuration variables to point to this. (there are defaults set
  for this)
* Site specific cluster parameters. We assume you are using SGE, you may
  wish to set:
    * *queue*: the queue to use on the cluster
    * *parallel_environment*: pe to use when submitting multicore jobs
    * *pe_queue*: queue to use for multicore jobs if different from above
    * *memory_resource*: comma separated list of resources to use when
      requesting memory for a job.
  Alternatively, run the pipeline with `--no-cluster` to run all jobs locally
  but be aware that this might take a very long time and a lot of memory.


Input
-----

The inputs are fastq.gz files. For the analysis of SRSF proteins the following
fastq files were downloaded from ENA:
    
+--------------+----------+
|*Sample*      |*Record*  |    
+--------------+----------+
|SRSF1-GFP-R1  |SRR2057564|
+--------------+----------+
|SRSF1-GFP-R2  |SRR2057565|
+--------------+----------+
|SRSF1-GFP-R3  |SRR2057566|
+--------------+----------+
|SRSF1-GFP-R4  |SRR2057567|
+--------------+----------+
|SRSF2-GFP-R1  |SRR2057568|
+--------------+----------+
|SRSF2-GFP-R2  |SRR2057569|
+--------------+----------+
|SRSF2-GFP-R3  |SRR2057570|
+--------------+----------+
|SRSF2-GFP-R4  |SRR2057571|
+--------------+----------+
|SRSF2-GFP-R5  |SRR2057572|
+--------------+----------+
|SRSF3-GFP-R1  |SRR2057573|
+--------------+----------+
|SRSF3-GFP-R2  |SRR2057574|
+--------------+----------+
|SRSF3-GFP-R3  |SRR2057575|
+--------------+----------+
|SRSF4-GFP-R1  |SRR2057576|
+--------------+----------+
|SRSF4-GFP-R2  |SRR2057577|
+--------------+----------+
|SRSF4-GFP-R3  |SRR2057578|
+--------------+----------+
|SRSF5-GFP-R1  |SRR2057579|
+--------------+----------+
|SRSF5-GFP-R2  |SRR2057580|
+--------------+----------+
|SRSF5-GFP-R3  |SRR2057581|
+--------------+----------+
|SRSF5-GFP-R4  |SRR2057582|
+--------------+----------+
|SRSF6-GFP-R1  |SRR2057583|
+--------------+----------+
|SRSF6-GFP-R2  |SRR2057584|
+--------------+----------+
|SRSF6-GFP-R3  |SRR2057585|
+--------------+----------+
|SRSF7-GFP-R1  |SRR2057586|
+--------------+----------+
|SRSF7-GFP-R2  |SRR2057587|
+--------------+----------+
|SRSF7-GFP-R3  |SRR2057588|
+--------------+----------+
|SRSF7-GFP-R4  |SRR2057589|
+--------------+----------+
|SRSF7-GFP-R5  |SRR2057590|
+--------------+----------+
|SRSF7-GFP-R6  |SRR2057591|
+--------------+----------+
|Nxf1-GFP-R1   |SRR2057592|
+--------------+----------+
|Nxf1-GFP-R2   |SRR2057593|
+--------------+----------+
|Nxf1-GFP-R3   |SRR2057594|
+--------------+----------+
|Control-GFP-R1|SRR2057595|
+--------------+----------+
|Control-GFP-R2|SRR2057596|
+--------------+----------+
|Control-GFP-R3|SRR2057597|
+--------------+----------+
|Control-GFP-R4|SRR2057598|
+--------------+----------+

These should be placed in the directory in which the pipeline is to be run. 
The pipeline expects these to be raw fastq files. That is that they contain the 
UMIs and the barcodes still on the 5' end  of the reads. This lite version of
iCLIP pipeline is expecting that all reads
associated with a sample are in the same fastq file.

In addition to the fastq files, a table of barcodes and samples is required as
sample_table.tsv. The data directory of the repository contains these for
the SRSF and TDP data. 

It has four columns:

The first contains the barcode including UMI bases, marked as Xs.
The second contains the barcode sequence without the UMI bases.
The third contains the sample name you'd like to use
The fourth contains the fastq files that contain reads from this sample

e.g.

NNNGGTTNN	GGTT	Control-GFP-R1	SRR12345678

Means that the sample FlipIn-FLAG-R1 should have reads in the fastq file
SRR12345678 is marked by the barcode GGTT and is embeded in the
UMI as NNNGGTTNN.

The file used to run the SRSF analysis is included in the data directory of the
repository and can be copied to the run directory.

Requirements
------------

The pipeline has been built to run in a linux environment. It may be possible to
run in OSX, but this is untested. The pipeline requires the following software 
to be installed (listed the versions used in this analysis):

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|CGAPipelines        | e6bb3be           |Pipelining infrastructure, mapping pipeline     |
|                    |                   |(http:/github.com/CGATOxford/CGATPipelines)     | 
+--------------------+-------------------+------------------------------------------------+
|CGAT                | 0.2.4             |Various                                         |
|                    |                   |(http:/github.com/CGATOxford/cgat)              |
+--------------------+-------------------+------------------------------------------------+
|Bowtie              | 1.1.2             |Mapping reads                                   |
+--------------------+-------------------+------------------------------------------------+
|FastQC              | 0.11.2            |Quality Control of demuxed reads                |
+--------------------+-------------------+------------------------------------------------+
|bedtools            | 2.22.0            |Interval manipulation                           |
+--------------------+-------------------+------------------------------------------------+
|samtools            | 1.3.1             |Read manipulation                               |
+--------------------+-------------------+------------------------------------------------+
|UMI-tools           | 0.0.2             |UMI manipulation                                |
+--------------------+-------------------+------------------------------------------------+
|reaper              | 13-100            |Used for demuxing and clipping reads            |
+--------------------+-------------------+------------------------------------------------+
|jupyter             | 4.1               |Running the statistical analysis and generating |
|                    |                   |figures                                         |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The for each deduping method there is a directory containing the deduped
BAM files, signifcant cluster calls  and statistics for each sample. e.g.
<adjacency.dir> contains the results for deduping using the adjacency method. 

The <plots> directory contains the image files for the figures used in the publication.

The <notebooks> directory contains html builds of the analysis notebooks exectuted
on the results of the analysis.


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import re

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping

import PipelineUMI

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.PARAMS["pipeline_src"] = os.path.dirname(__file__)
P.getParameters(
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
if '%' in PARAMS["sample_table"]:
    PARAMS["sample_table"] = PARAMS["sample_table"] % PARAMS["pipeline_src"]
if '%' in PARAMS["annotations_gtf"]:
    PARAMS["annotations_gtf"] = PARAMS["annotations_gtf"] % PARAMS["pipeline_src"]
if '%' in PARAMS["annotations_contigs"]:
    PARAMS["annotations_contigs"] = PARAMS["annotations_contigs"] % PARAMS["pipeline_src"]

###################################################################
###################################################################
###################################################################
## WORKER TASKS
###################################################################
# Read preparation
###################################################################
@jobs_limit(1, "db")
@transform(PARAMS["sample_table"], suffix(".tsv"), ".load")
def loadSampleInfo(infile, outfile):

    P.load(infile, outfile,
           options="--header-names=format,barcode,track,lanes -i barcode -i track")


###################################################################
@follows(mkdir("demux_fq"))
@transform("*.fastq.gz", regex("(.+).fastq.gz"),
           r"demux_fq/\1.fastq.umi_trimmed.gz")
def extractUMI(infile, outfile):
    ''' Remove UMI from the start of each read and add to the read
    name to allow later deconvolving of PCR duplicates '''

    statement=''' zcat %(infile)s
                | umi_tools extract
                        --bc-pattern=NNNXXXXNN
                        -L %(outfile)s.log
                | gzip > %(outfile)s '''

    P.run()


###################################################################
@jobs_limit(1, "db")
@transform(extractUMI, suffix(".fastq.umi_trimmed.gz"),
           "umi_stats.load")
def loadUMIStats(infile, outfile):
    ''' load stats on UMI usage from the extract_umi log into the
    database '''

    infile = infile + ".log"
    P.load(infile, outfile, "-i sample -i barcode -i UMI")


###################################################################
@transform("*.fastq.gz",
           regex("(.+).fastq.gz"),
           add_inputs(PARAMS["sample_table"]),
           r"\1_reaper_metadata.tsv")
def generateReaperMetaData(infile, outfile):
    '''Take the sample_table and use it to generate a metadata table
    for reaper in the correct format '''

    adaptor_5prime = "AGATCGGAAGAGCGACGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    adaptor_3prime = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"

    outlines = []
    lane = P.snip(infile[0], ".fastq.gz")
    for line in IOTools.openFile(infile[1]):
        fields = line.split("\t")
        barcode = fields[1]
        lanes = fields[-1].strip().split(",")
        if lane in lanes:
            outlines.append([barcode, adaptor_3prime, adaptor_5prime, "-"])

    header = ["barcode", "3p-ad", "tabu", "5p-si"]
    IOTools.writeLines(outfile, outlines, header)


###################################################################
@follows(loadUMIStats, generateReaperMetaData)
@subdivide(extractUMI, regex(".+/(.+).fastq.umi_trimmed.gz"),
           add_inputs(r"\1_reaper_metadata.tsv", PARAMS["sample_table"]),
           r"demux_fq/*_\1.fastq.gz")
def demux_fastq(infiles, outfiles):
    '''Demultiplex each fastq file into a seperate file for each
    barcode/UMI combination'''

    infile, meta, samples = infiles
    track = re.match(".+/(.+).fastq.umi_trimmed.gz", infile).groups()[0]

    statement = '''reaper -geom 5p-bc
                          -meta %(meta)s
                          -i <( zcat %(infile)s | sed 's/ /_/g')
                          --noqc
                          -3p-head-to-tail 2 
                          -3p-prefix 6/2/1
                          -basename demux_fq/%(track)s_
                          -clean-length 15 > %(track)s_reapear.log;
                   checkpoint;
                   rename _. _ demux_fq/*clean.gz;
                 '''

    for line in IOTools.openFile(samples):
        line = line.split("\t")
        bc, name, lanes = line[1:]
        name = name.strip()
           
        if track in lanes.strip().split(","):
            statement += '''checkpoint;
                            mv demux_fq/%(track)s_%(bc)s.clean.gz
                               demux_fq/%(name)s_%(track)s.fastq.gz; ''' % locals()

    P.run()


###################################################################
@follows(mkdir("fastqc"))
@transform(demux_fastq, regex(".+/(.+).fastq(.*)\.gz"),
           r"fastqc/\1\2.fastqc")
def qcDemuxedReads(infile, outfile):
    ''' Run fastqc on the post demuxing and trimmed reads'''

    m = PipelineMapping.FastQc(nogroup=False, outdir="fastqc")
    statement = m.build((infile, ), outfile)
    exportdir = "fastqc"
    P.run()


###################################################################
@follows(demux_fastq, qcDemuxedReads, loadUMIStats)
def PrepareReads():
    pass


###################################################################
# Mapping
###################################################################
@follows(mkdir("mapping.dir"), demux_fastq)
@transform(demux_fastq,
           regex(".+/(.+)_(.+).fastq.gz"),
           r"mapping.dir/\1.bam")
def run_mapping(infile, outfile):
    ''' Map reads using the selected read mapper '''

    job_threads = PARAMS["bowtie_threads"]
    job_memory = PARAMS["bowtie_memory"]

    m = PipelineMapping.Bowtie(
           executable="bowtie",
           tool_options=PARAMS["bowtie_options"],
           strip_sequence=0)

    genome = PARAMS["bowtie_genome"]
    reffile = os.path.join(PARAMS["bowtie_index_dir"],
                           PARAMS["bowtie_genome"] + ".fa")

    statement = m.build((infile,), outfile)

    P.run()


###################################################################
# Deduping, Counting, etc
###################################################################
# dedup methods
METHODS = ["unique", "cluster", "percentile", "adjacency", "directional"]


###################################################################
@originate(["dedup_%s.sentinal" % method for method in METHODS])
def dedup_method_sentinals(outfile):
    ''' make sentinal files for each dedup method so can be
    called with @product, and create output directories'''

    os.mkdir(P.snip(outfile, ".sentinal") + ".dir")
    P.touch(outfile)


###################################################################
@product(run_mapping,
         formatter(".+/(?P<TRACK>.+).bam"),
         dedup_method_sentinals,
         formatter("dedup_(?P<method>.+)\.sentinal"),
         "dedup_{method[1][0]}.dir/{TRACK[0][0]}.bam",
         ["{basename[0][0]}", "{method[1][0]}"])
def dedup_bams(infiles, outfile, params):
    '''run umi_tools dedup on all tracks using all methods'''

    track, method = params

    infile, outdir = infiles
    outdir = P.snip(outdir, ".sentinal") + ".dir"
    track = os.path.join(outdir, track)

    job_memory = "21G"

    if method == "cluster":
        further_stats = "--further-stats"
    else:
        further_stats = ""

    statement = []
    statement.append('''umi_tools dedup
                         --method=%(method)s
                         --output-stats=%(track)s
                         %(further_stats)s
                          -I %(infile)s
                          -S %(track)s.unsorted.bam
                          -L %(track)s.log ''')

    statement.append("samtools sort %(track)s.unsorted.bam -O bam -T %(track)s > %(track)s.bam")
    statement.append("samtools index %(track)s.bam")
    statement.append("rm %(track)s.unsorted.bam")
    statement = "; checkpoint;".join(statement)
    P.run()


###################################################################
@jobs_limit(1, "db")
@transform(dedup_bams,
           regex("(.+).dir/(.+).bam"),
           inputs(r"\1.dir/\2_edit_distance.tsv"),
           r"\1.dir/\2_\1_edit_distance.load")
def loadEditDistances(infile, outfile):
    '''Load distribtuions of edit distances as output by umi_tools dedup'''
    load_smt = P.build_load_statement(
        P.toTable(outfile), options="-i edit_distance")
    statement = ''' sed s/unique/_unique/g %(infile)s
                 | %(load_smt)s > %(outfile)s '''
    P.run()


###################################################################
@jobs_limit(1, "db")
@collate(dedup_bams,
         formatter("dedup_cluster.dir/(?P<Track>.+).bam"),
         inputs(r"dedup_cluster.dir/{Track[0]}_topologies.tsv"),
         r"topologies.load")
def load_topologies(infiles, outfile):
    '''Load the topologies distribution - only output if method was 
    cluster as it will be the same irrespective of the network
    method used'''

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_topologies.tsv",
                         has_titles=False,
                         header="track,category,count")


###################################################################
@jobs_limit(1, "db")
@collate(dedup_bams,
         formatter("dedup_cluster.dir/(?P<Track>.+).bam"),
         inputs(r"dedup_cluster.dir/{Track[0]}_nodes.tsv"),
         r"node_counts.load")
def load_node_counts(infiles, outfile):
    '''Load the number of counts per cluster distribution - only
    output if method was cluster as it will be the same irrespective
    of the network method used'''

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_nodes.tsv",
                         has_titles=False,
                         header="track,category,count")


##################################################################
@collate(dedup_bams,
         regex("(.+/.+-.+)-(R[0-9]+)(.*).bam"),
         r"\1-agg\3.bam")
def mergeBamsByRep(infiles, outfile):
    '''Merge all replicates of the same experiment together'''

    statement = '''samtools merge -f %(outfile)s %(infiles)s;
                   checkpoint;
                   samtools index %(outfile)s'''
    if len(infiles) == 1:
        IOTools.cloneFile(infiles[0], outfile)
        IOTools.cloneFile(infiles[0]+".bai", outfile+".bai")
    else:
        infiles = " ".join(infiles)
        P.run()


###################################################################
@merge([mergeBamsByRep, dedup_bams, "mapping.dir/*.bam"],
       "read_counts.tsv")
def count_bams(infiles, outfile):
    '''Count the number of alignments both pre and post dedup'''

    outlines = []

    for infile in infiles:

        method = re.match("dedup_(.+).dir\/.+", infile)
        if method:
            method = method.groups()[0]
        else:
            method = "none"
            
        track = re.search("([^/]+).bam", infile).groups()[0]

        statement = '''samtools idxstats %(infile)s
                 | awk '{sum+=$3} END{print sum}' '''

        count, _ = P.execute(statement)

        outlines.append([method, track, count.strip()])

    IOTools.writeLines(outfile,
                       outlines, header=["method", "track", "count"])


###################################################################
@jobs_limit(1, "db")
@transform(count_bams, suffix(".tsv"), ".load")
def load_read_counts(infile, outfile):

    P.load(infile, outfile, options="-i method -i track")


###################################################################
@follows(loadEditDistances,
         load_topologies,
         load_node_counts,
         load_read_counts)
def get_dedup_stats():
    pass


###################################################################
@follows(mkdir("plots"), mkdir("notebooks"),
         get_dedup_stats)
@transform([os.path.join(PARAMS["pipeline_src"], x)
            for x in ["notebooks/Fig1b_2a-Edit_distance_dists.ipynb",
                      "notebooks/Fig1c-Read_depth_and_errors.ipynb"]],
           regex(".+/(notebooks/.+).ipynb"),
           r"\1.html")
def runNotebooks1(infile, outfile):
    '''Run notebooks and generate figures upto this point'''

    ipynb_file = os.path.basename(infile)

    statement = '''cp %(infile)s %(ipynb_file)s; checkpoint;
                   ipython nbconvert --to=html
                                     --output=%(outfile)s
                                     --execute
                                     %(ipynb_file)s > %(outfile)s.log'''

    P.run()


###################################################################
# Analyses and clusters etc
###################################################################
@transform([dedup_bams, mergeBamsByRep],
           regex("(.+/.+).bam"),
           add_inputs(PARAMS["annotations_gtf"]),
           r"\1.clusters.bedgraph.gz")
def call_clusters_by_rand(infiles, outfile):
    '''Use randomisation within a gene to call significantly crosslinked
    bases'''

    bamfile, annotation = infiles

    job_threads = 6
    job_memory = "0.5G"
    statement = '''python %(pipeline_src)s/iCLIPlib/scripts/significant_bases_by_randomisation.py
                      -b %(bamfile)s
                      -I %(annotation)s
                      -p 6
                      -S %(outfile)s
                      -L %(outfile)s.log'''

    P.run()


###################################################################
@jobs_limit(1, "db")
@merge(call_clusters_by_rand,
       "cluster_counts.load")
def loadClusterCounts(infiles, outfile):
    '''Find the number of signficant clusters found in each sample'''

    tmp = P.getTempFilename(shared=True)
    results = []
    for infile in infiles:
        count = IOTools.getNumLines(infile)
        method, track = re.match(
            "dedup_(.+).dir/(.+)\.clusters.bedgraph", infile).groups()
        results.append((method, track, count))
        
    IOTools.writeLines(tmp, results, header=["method", "track", "count"])

    P.load(tmp, outfile)
    os.unlink(tmp)


###################################################################
@transform(call_clusters_by_rand,
           regex("(.+).clusters.bedgraph.gz"),
           add_inputs(r"\1.bam"),
           r"\1.sig_bases.bedgraph.gz")
def get_sig_bases(infiles, outfile):
    '''significant_bases_by_randomiasation returns pvalues
    but for the next step we need counts. retrieve the significant
    bases and get their tag counts'''

    sig_file, bamfile = infiles
    PipelineUMI.getSigHeights(sig_file,
                              bamfile,
                              outfile,
                              submit=True)


###################################################################
@transform(call_clusters_by_rand,
           regex("(.+/.+).clusters.bedgraph.gz"),
           r"\1.merged_clusters.bed.gz")
def merge_adjacent_clusters(infile, outfile):
    '''Merge bases called as significant if their territories overlap'''

    statement = '''bedtools slop -b 15 -i %(infile)s -g %(annotations_contigs)s
                 | sort -k1,1 -k2,2n
                 | bedtools merge -i -
                 | gzip > %(outfile)s'''

    P.run()


###################################################################
@follows(loadClusterCounts, get_sig_bases, merge_adjacent_clusters)
def clusters():
    pass


###################################################################
# Analysis of correlation of Exons
###################################################################
@transform(PARAMS["annotations_gtf"],
           regex(".+gtf.gz"),
           "intersected_exons.gtf.gz")
def intersect_exons(infile, outfile):
    '''Take each gene and where exons from different transcripts
    overlap, take the intersection of those exons'''
    statement = ''' python %(scriptsdir)s/gtf2gtf.py
                           -I %(infile)s
                           --method=intersect-transcripts
                           
                           -S %(outfile)s '''
    P.run()


###################################################################
@transform([dedup_bams, mergeBamsByRep],
           suffix(".bam"),
           add_inputs(intersect_exons),
           ".exon_count.tsv.gz")
def count_exons(infiles, outfile):
    '''Count the number of clip tags in each exon'''
    bamfile, gtffile = infiles

    statement = '''python %(pipeline_src)s/iCLIPlib/scripts/count_clip_sites.py
                          -I %(gtffile)s
                          %(bamfile)s
                           -f exon
                          -S %(outfile)s '''

    P.run()


###################################################################
@jobs_limit(1, "db")
@merge(count_exons, "exon_counts.load")
def load_exon_counts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="dedup_(.+).dir/(.+).exon_count.tsv.gz",
                         cat="method,track",
                         options="-i method -i track -i gene_id")


###################################################################
@transform(get_sig_bases,
           suffix("sig_bases.bedgraph.gz"),
           add_inputs(intersect_exons),
           "exon_sig_count.tsv.gz")
def count_sig_bases_over_exons(infiles, outfile):
    '''Count the number of tags in significant bases that
    overlap with each exon'''

    bedgraph, bed = infiles

    statement = '''
                   python %(scriptsdir)s/gff2bed.py --is-gtf -I %(bed)s -L %(outfile)s.log
                 | sort -k1,1 -k2,2n
                 | bedtools map -a - -b %(bedgraph)s -c 4 -o sum -null 0
                 | cut -f4,7
                 | awk 'BEGIN{OFS="\\t"; i=0} {$1=$1 "_" (i++); print $0}' 
                 | awk '$2 > 0'
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@jobs_limit(1, "db")
@merge(count_sig_bases_over_exons, "sig_exon_counts.load")
def load_sig_exon_counts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="dedup_(.+).dir/(.+).exon_sig_count.tsv.gz",
                         cat="method,track",
                         has_titles=False,
                         header="method,track,name,count",
                         options="-i method -i track")


###################################################################
@follows(load_sig_exon_counts, load_exon_counts)
def exon_level_correlation():
    pass


###################################################################
# Base level reproducibility
###################################################################
@collate(dedup_bams,
         regex("(.+/.+)-(R.).bam"),
         r"\1.reproducibility.tsv")
def calculate_base_level_reproducibility(infiles, outfile):
    '''Find the number of bases in one sample clipping in others'''

    infiles = " ".join([infile for infile in infiles if re.search("R[123]", infile)])
    statement = '''python %(pipeline_src)s/iCLIPlib/scripts/calculateiCLIPReproducibility.py
                       %(infiles)s
                       -m 2
                       -S %(outfile)s
                       -L %(outfile)s.log'''
    P.run()


###################################################################
@jobs_limit(1, "db")
@merge(calculate_base_level_reproducibility, "base_level_reproducibility.load")
def load_base_level_reproducibility(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="dedup_(.+).dir/.+.rep",
                         cat="method",
                         options="-i method")


###################################################################
@follows(load_base_level_reproducibility)
def base_level_reproducibility():
    pass


###################################################################
@follows(mkdir("plots"), mkdir("notebooks"),
         get_dedup_stats,
         base_level_reproducibility,
         exon_level_correlation,
         clusters)
@transform([os.path.join(PARAMS["pipeline_src"], "notebooks/Fig*.ipynb"),
            os.path.join(PARAMS["pipeline_src"], "notebooks/Examine*.ipynb")],
           regex(".+/(notebooks/.+).ipynb"),
           r"\1.html")
def runNotebooks2(infile, outfile):
    '''Run all notebooks and generate figures'''

    ipynb_file = os.path.basename(infile)

    statement = '''cp %(infile)s %(ipynb_file)s; checkpoint;
                   ipython nbconvert --to=html
                                     --output=%(outfile)s
                                     --ExecutePreprocessor.timeout=-1
                                     --execute
                                     %(ipynb_file)s > %(outfile)s.log'''

    P.run()


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, run_mapping,
         get_dedup_stats,
         load_read_counts,
         clusters,
         exon_level_correlation,
         base_level_reproducibility,
         runNotebooks2)
def full():
    pass


@follows( mkdir( "report" ))
def build_report():
    '''build report from scratch.'''

    try:
        os.symlink(os.path.abspath("conf.py"),
                   os.path.join(
                       os.path.abspath("mapping.dir"), "conf.py"))
    except OSError as e:
        E.warning(str(e))

    E.info("Running mapping report build from scratch")
#    statement = '''cd mapping.dir;
#                   python %(scripts_dir)s/CGATPipelines/pipeline_mapping.py
#                   -v5 -p1 make build_report '''
#    P.run()
    E.info("starting report build process from scratch")
    P.run_report(clean = True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
