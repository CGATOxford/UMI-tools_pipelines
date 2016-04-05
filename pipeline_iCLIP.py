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
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use
CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini`
file  (see :ref:`PipelineDocumenation`). To start with, use the files supplied
with the :ref:`Example` data.

Input
-----

The inputs should be a fastq file. The pipeline expects these to be raw fastq
files. That is that they contain the UMIs and the barcodes still on the 5' end
of the reads. It is expected that these will not be demultiplexed, although
demultiplex fastqs will also work.

In addition to the fastq files, a table of barcodes and samples is required as
sample_table.tsv.

It has four columns:

The first contains the barcode including UMI bases, marked as Xs.
The second contains the barcode sequence without the UMI bases.
The third contains the sample name you'd like to use
The fourth contains the fastq files that contain reads from this sample

e.g.

NNNGGTTNN	GGTT	FlipIn-FLAG-R1	hiseq7,hiseq8,miseq

Means that the sample FlipIn-FLAG-R1 should have reads in the fastq files
hiseq7, hiseq8 and miseq, is marked by the barcode GGTT and is embeded in the
UMI as NNNGGTTNN.

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|CGAPipelines        |                   |Pipelining infrastructure, mapping pipeline     |
+--------------------+-------------------+------------------------------------------------+
|CGAT                | >=0.2.4           |Various                                         |
+--------------------+-------------------+------------------------------------------------+
|Bowtie              | >=1.1.2           |Filtering out PhiX reads                        |
+--------------------+-------------------+------------------------------------------------+
|FastQC              | >=0.11.2          |Quality Control of demuxed reads                |
+--------------------+-------------------+------------------------------------------------+
|STAR or HiSat       |                   |Spliced mapping of reads                        |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |Interval manipulation                           |
+--------------------+-------------------+------------------------------------------------+
|samtools            |                   |Read manipulation                               |
+--------------------+-------------------+------------------------------------------------+
|subread             |                   |FeatureCounts read quantifiacation              |
+--------------------+-------------------+------------------------------------------------+
|MEME                |>=4.9.1            |Motif finding                                   |
+--------------------+-------------------+------------------------------------------------+
|DREME               |>=4.9.1            |Motif finding                                   |
+--------------------+-------------------+------------------------------------------------+
|bedGraphToBigWig    |                   |Converstion of results to BigWig                |
+--------------------+-------------------+------------------------------------------------+
|reaper              | 13-100            |Used for demuxing and clipping reads            |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

As well as the report, clusters, as BED files are in the clusters.dir directory,
and traces as bigWig files are in the bigwig directory. Both of these are exported
as a UCSU genome browser track hub in the export directory. 

Example
=======

Example data and configuration is avaiable in example_data.tar.gz


Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys
import glob
import os
import re

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(PARAMS["annotations_dir"],
                                      "pipeline_annotations.py")

###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################
@transform("sample_table.tsv", suffix(".tsv"), ".load")
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
                        --bc-pattern=%(reads_bc_pattern)s
                        -L %(outfile)s.log
                | gzip > %(outfile)s '''

    P.run()


###################################################################
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
           add_inputs("sample_table.tsv"),
           r"\1_reaper_metadata.tsv")
def generateReaperMetaData(infile, outfile):
    '''Take the sample_table and use it to generate a metadata table
    for reaper in the correct format '''

    adaptor_5prime = PARAMS["reads_5prime_adapt"]
    adaptor_3prime = PARAMS["reads_3prime_adapt"]

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
           add_inputs(r"\1_reaper_metadata.tsv", "sample_table.tsv"),
           r"demux_fq/*.fastq.gz")
def demux_fastq(infiles, outfiles):
    '''Demultiplex each fastq file into a seperate file for each
    barcode/UMI combination'''

    infile, meta, samples = infiles
    track = re.match(".+/(.+).fastq.umi_trimmed.gz", infile).groups()[0]

    statement = '''reaper -geom 5p-bc
                          -meta %(meta)s
                          -i <( zcat %(infile)s | sed 's/ /_/g')
                          --noqc
                          %(reads_reaper_options)s
                          -basename demux_fq/%(track)s_
                          -clean-length %(reads_min_length)s > %(track)s_reapear.log;
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
                            demux_fq/%(name)s.fastq.gz; ''' % locals()

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
           regex(".+/(.+).fastq.gz"),
           r"mapping.dir/\1.bam")
def run_mapping(infile, outfiles):
    ''' run the mapping target of the mapping pipeline '''

    if PARAMS["mapper"] == "bowtie":
        job_threads = PARAMS["bowtie_threads"]
        job_memory = PARAMS["bowtie_memory"]

        m = PipelineMapping.Bowtie(
            executable=PARAMS["bowtie_executable"],
            tool_options=PARAMS["bowtie_options"],
            strip_sequence=PARAMS["strip_sequence"])
        
        reffile = os.path.join(PARAMS["bowtie_index_dir"],
                               PARAMS["genome"] + ".fa")
        statement = m.build((infile,), outfile)

    elif["mapper"] == "star":
        job_threads = PARAMS["star_threads"]
        job_memory = PARAMS["star_memory"]

        star_mapping_genome = PARAMS["star_genome"] or PARAMS["genome"]

        m = PipelineMapping.STAR(
            executable=P.substituteParameters(**locals())["star_executable"],
            strip_sequence=PARAMS["strip_sequence"])

        statement = m.build((infile,), outfile)
    else:
        raise ValueError("Mapper '%s' not implemented" % PARAMS["mapper"])

    P.run()


###################################################################
# Deduping, Counting, etc
###################################################################
@follows(mkdir("deduped.dir"), run_mapping)
@transform(run_mapping, regex("(.+)/(.+)\.([^\.]+)\.bam"),
           r"deduped.dir/\2.bam")
def dedup_alignments(infile, outfile):
    ''' Deduplicate reads, taking UMIs into account'''

    outfile = P.snip(outfile, ".bam")

    job_memory = "7G"
    statement = ''' umi_tools dedup
                    %(dedup_options)s
                    --output-stats=%(outfile)s
                    -I %(infile)s
                    -S %(outfile)s.tmp.bam
                    -L %(outfile)s.log;
 
                    checkpoint;

                    samtools sort %(outfile)s.tmp.bam %(outfile)s;
                   
                    checkpoint;

                    samtools index %(outfile)s.bam ;
 
                    checkpoint;

                    rm %(outfile)s.tmp.bam'''

    P.run()


###################################################################
@transform(dedup_alignments, suffix(".bam"), ".umi_stats.tsv.gz")
def deduped_umi_stats(infile, outfile):
    ''' calculate histograms of umi frequencies '''

    statement = '''python %(project_src)s/umi_hist.py
                           -I %(infile)s
                           -L %(outfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
@merge(deduped_umi_stats, "dedup_umi_stats.load")
def loadDedupedUMIStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).umi_stats.tsv.gz",
                         cat="track",
                         options="-i track -i UMI")


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, mapping)
def full():
    pass


@follows( mkdir( "report" ), createViewMapping)
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


@follows(mkdir("report"), createViewMapping)
def update_report():
    '''update report.'''

    E.info("Updating Mapping Reported")
    statement = '''cd mapping.dir;
                   python %(pipelinedir)s/pipeline_mapping.py
                   -v5 -p1 make update_report '''
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
