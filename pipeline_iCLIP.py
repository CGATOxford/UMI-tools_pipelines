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
#################################################################################
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

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMotifs as PipelineMotifs
import CGATPipelines.PipelineRnaseq as PipelineRnaseq
import PipelineiCLIP
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
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

PipelineMotifs.PARAMS = PARAMS
PipelineiCLIP.PARAMS = PARAMS
PipelineiCLIP.PARAMS_ANNOTATIONS = PARAMS_ANNOTATIONS

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3)
for line in IOTools.openFile("sample_table.tsv"):
    track = line.split("\t")[2]
    TRACKS.tracks.append(PipelineTracks.Sample3(filename=track))



###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' \
                % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################
@transform("*.fastq.gz", regex("(.+).fastq.gz"),
           add_inputs(os.path.join(PARAMS["bowtie_index_dir"],
                                   PARAMS["phix_genome"]+".fa")),
           r"\1.fastq.clean.gz")
def filterPhiX(infiles, outfile):
    ''' Use mapping to bowtie to remove any phiX mapping reads '''

    infile, reffile = infiles
    outfile = P.snip(outfile, ".gz")
    bam_out = P.snip(infile, ".fastq.gz") + ".phix.bam"

    job_threads = PARAMS["phix_bowtie_threads"]
    job_memory = PARAMS["phix_bowtie_memory"]
    options = PARAMS["phix_bowtie_options"] + " --un %s" % outfile
    genome = PARAMS["phix_genome"]
    bowtie_threads = PARAMS["phix_bowtie_threads"]

    m = PipelineMapping.Bowtie(executable=PARAMS["phix_bowtie_exe"],
                               strip_sequence=False,
                               remove_non_unique=False,
                               tool_options=options)

    statement = m.build((infile,), bam_out)
    statement += "checkpoint; gzip %(outfile)s"

    P.run()


@transform("sample_table.tsv", suffix(".tsv"), ".load")
def loadSampleInfo(infile, outfile):

    P.load(infile, outfile,
           options="--header-names=format,barcode,track,lanes -i barcode -i track")
###################################################################
@follows(mkdir("demux_fq"))
@transform(filterPhiX, regex("(.+).fastq.clean.gz"),
           r"demux_fq/\1.fastq.umi_trimmed.gz")
def extractUMI(infile, outfile):
    ''' Remove UMI from the start of each read and add to the read
    name to allow later deconvolving of PCR duplicates '''

    statement=''' zcat %(infile)s
                | python %(project_src)s/UMI-tools/extract_umi.py
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
@transform(filterPhiX,
           regex("(.+).fastq.clean.gz"),
           add_inputs("sample_table.tsv"),
           r"\1_reaper_metadata.tsv")
def generateReaperMetaData(infile, outfile):
    '''Take the sample_table and use it to generate a metadata table
    for reaper in the correct format '''

    adaptor_5prime = PARAMS["reads_5prime_adapt"]
    adaptor_3prime = PARAMS["reads_3prime_adapt"]

    outlines = []
    lane = P.snip(infile[0], ".fastq.clean.gz")
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
       r"demux_fq/*_\1.fastq*gz")
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
        if PARAMS["reads_paired"]:
            ext = "fastq.1.gz"
        else:
            ext = "fastq.gz"
        if track in lanes.strip().split(","):
            statement += '''checkpoint;
                         mv demux_fq/%(track)s_%(bc)s.clean.gz
                            demux_fq/%(name)s_%(track)s.%(ext)s; ''' % locals()

    P.run()


###################################################################
@active_if(PARAMS["reads_paired"]==1)
@transform("*.fastq.2.gz", suffix(".fastq.2.gz"),
           ".fastq.reaped.2.gz")
def reapRead2(infile,outfile):

    track = P.snip(outfile,".fastq.reaped.2.gz")
    statement = ''' reaper -geom no-bc
                           -3pa %(reads_3prime_adapt)s
                           -i %(infile)s
                           -basename %(track)s
                           --noqc
                           -clean-length 15 > %(track)s_pair_reaper.log;
                    checkpoint;
                    mv %(track)s.lane.clean.gz %(outfile)s '''

    P.run()


###################################################################
@active_if(PARAMS["reads_paired"]==1)
@follows(mkdir("reconciled.dir"), reapRead2)
@transform(demux_fastq,
           regex(".+/(.+)_(.+).fastq.1.gz"),
           add_inputs(r"\2.fastq.2.gz"),
           [r"reconciled.dir/\1_\2.fastq.2.gz",
            r"reconciled.dir/\1_\2.fastq.1.gz"])
def reconsilePairs(infiles, outfiles):
    ''' Pull reads read 2 file that are present in each read 1 file '''

    track = P.snip(os.path.basename(infiles[0]), ".fastq.1.gz")
    infiles = " ".join(infiles)
    job_options = "-l mem_free=1G"
    statement = '''python %(scripts_dir)s/fastqs2fastqs.py
                          --method=reconcile
                          --id-pattern-1='(.+)_.+_[ATGC]+'
                          --output-filename-pattern=reconciled.dir/%(track)s.fastq.%%s.gz
                           %(infiles)s > reconciled.dir/%(track)s.log '''

    P.run()


###################################################################
@follows(mkdir("fastqc"))
@transform(demux_fastq, regex(".+/(.+).fastq(.*)\.gz"),
           r"fastqc/\1\2.fastqc")
def qcDemuxedReads(infile, outfile):
    ''' Run fastqc on the post demuxing and trimmed reads'''

    m = PipelineMapping.FastQc(nogroup=False, outdir="fastqc")
    statement = m.build((infile,),outfile)
    exportdir = "fastqc"
    P.run()


###################################################################
@transform(qcDemuxedReads, regex("(.+)/(.+)\.fastqc"),
           inputs(r"\1/\2_fastqc/fastqc_data.txt"),
           r"\1/\2_length_distribution.tsv")
def getLengthDistribution(infile, outfile):
    ''' Parse length distribution out of the fastqc results '''

    statement = '''
       sed -e '/>>Sequence Length Distribution/,/>>END_MODULE/!d' %(infile)s
     | grep -P '^[0-9]+'
     | sed '1istart\\tend\\tcount'
     | sed 's/-/\\t/' > %(outfile)s '''

    P.run()


###################################################################
@merge(getLengthDistribution, "read_length_distribution.load")
def loadLengthDistribution(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_length_distribution.tsv",
                         options="-i start -i end")


###################################################################
@follows(demux_fastq,qcDemuxedReads, loadUMIStats, loadSampleInfo,
         loadLengthDistribution)
def PrepareReads():
    pass


###################################################################
# Mapping
###################################################################
def mapping_files():

    infiles = glob.glob("%s/*.fastq*gz" % PARAMS["input"])
    infiles = [infile for infile in infiles if ".fastq.umi_trimmed.gz" not in infile]
    outfiles = set(
        ["mapping.dir/%(mapper)s.dir/%(track)s.%(mapper)s.bam"
         % {"mapper": PARAMS["mappers"],
            "track": re.match("%s/(.+).fastq.*gz"
                              % PARAMS["input"], infile).groups()[0]}
         for infile in infiles
         if "umi_trimmed" not in infile])

    yield (infiles, outfiles)


###################################################################
@follows(mkdir("mapping.dir"), demux_fastq)
@files(mapping_files)
def run_mapping(infiles, outfiles):
    ''' run the mapping target of the mapping pipeline '''

    to_cluster = False
    statement = ''' ln -f pipeline.ini mapping.dir/pipeline.ini;
                    checkpoint;
                    cd mapping.dir;
                    nice python %(pipelinedir)s/pipeline_mapping.py
                    make mapping
                    -v5 -p%(pipeline_mapping_jobs)s  '''

    P.run()


###################################################################
@follows(run_mapping)
@collate("mapping.dir/*.dir/*-*-*_*.bam",
         regex("(.+)/([^_]+\-.+\-[^_]+)_(.+)\.([^.]+)\.bam"),
         r"\1/merged_\2.\4.bam")
def mergeBAMFiles(infiles, outfile):
    '''Merge reads from the same library, run on different lanes '''

    if len(infiles) == 1:
        P.clone(infiles[0], outfile)
        P.clone(infiles[0]+".bai", outfile+".bai")
        return

    infiles = " ".join(infiles)
    statement = ''' samtools merge %(outfile)s %(infiles)s >& %(outfile)s.log'''
    P.run()


###################################################################
@transform(mergeBAMFiles, suffix(".bam"), ".bam.bai")
def indexMergedBAMs(infile, outfile):
    ''' Index the newly merged BAM files '''

    statement = ''' samtools index %(infile)s '''
    P.run()


###################################################################
@merge(indexMergedBAMs, "mapping.sentinal")
def mapping_qc(infiles, outfile):
    ''' run mapping pipeline qc targets '''

    to_cluster = False

    statement = '''cd mapping.dir;
                   nice python %(pipelinedir)s/pipeline_mapping.py
                   make qc -v5 -p%(pipeline_mapping_jobs)s '''
    P.run()

    P.touch("mapping.sentinal")


###################################################################
@follows(mapping_qc)
@originate("mapping.dir/geneset.dir/reference.gtf.gz")
def buildReferenceGeneSet(outfile):
    
    to_cluster = False

    statement = '''cd mapping.dir;
                   nice python %(pipelinedir)s/pipeline_mapping.py
                        make buildReferenceGeneSet
                        -v5 -p%(pipeline_mapping_jobs)s '''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.context.bed.gz")
def generateContextBed(infile, outfile):
    ''' Generate full length primary transcript annotations to count
    mapping contexts '''

    genome = os.path.join(PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_contigs"])
    statement = ''' zcat %(infile)s
                  | awk '$3=="exon"'
                  | python %(scripts_dir)s/gtf2gtf.py
                    --method=exons2introns
                    
                     -L %(outfile)s.log
                  | awk 'BEGIN{FS="\\t";OFS="\\t"} {$2="intron"; print}'
                  | gzip > %(outfile)s.tmp.gtf.gz;

                  checkpoint;

                  zcat %(infile)s %(outfile)s.tmp.gtf.gz
                  | awk '$3=="exon" || $3=="intron"'
                  | python %(scripts_dir)s/gff2bed.py
                    --set-name=source
                     -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
                    > %(outfile)s.tmp.bed;
                    
                    checkpoint;
                    
                    bedtools complement -i %(outfile)s.tmp.bed -g %(genome)s
                  | awk 'BEGIN{OFS="\\t"} {print ($0 "\\tnone\\t0\\t.")}'
                    > %(outfile)s.tmp2.bed;

                    checkpoint;

                    cat %(outfile)s.tmp.bed %(outfile)s.tmp2.bed
                  | sort -k1,1 -k2,2n
                  | gzip > %(outfile)s;
                  
                    checkpoint;
                  
                    rm %(outfile)s.tmp.bed %(outfile)s.tmp2.bed %(outfile)s.tmp.gtf.gz '''

    P.run()


###################################################################
@transform(generateContextBed, suffix(".context.bed.gz"),
           ".context_interval_stats.tsv.gz")
def getContextIntervalStats(infile, outfile):
    ''' Generate length stastics on context file '''

    statement = ''' python %(scripts_dir)s/bed2stats.py
                            --aggregate-by=name
                            -I %(infile)s
                    | gzip > %(outfile)s '''

    P.run()



###################################################################
@transform(getContextIntervalStats, suffix(".tsv.gz"),
           ".load")
def loadContextIntervalStats(infile, outfile):

    P.load(infile, outfile)


###################################################################
@transform(mapping_qc, regex("(.+)"),
           "mapping.dir/view_mapping.load")
def createViewMapping(infile, outfile):
    ''' Create tables neccessary for mapping report '''

    to_cluster = False
    statement = '''cd mapping.dir;
                   nice python %(scripts_dir)s/../../CGATPipelines/CGATPipelines/pipeline_mapping.py
                   make createViewMapping -v5 -p1 '''
    P.run()


###################################################################
@follows(mapping_qc,loadContextIntervalStats )
def mapping():
    pass


###################################################################
# Deduping, Counting, etc
###################################################################
@follows(mkdir("deduped.dir"), run_mapping)
@transform(indexMergedBAMs, regex("(.+)/merged_(.+)\.([^\.]+)\.bam.bai"),
           inputs(r"\1/merged_\2.\3.bam"),
           r"deduped.dir/\2.bam")
def dedup_alignments(infile, outfile):
    ''' Deduplicate reads, taking UMIs into account'''

    outfile = P.snip(outfile, ".bam")

    job_memory="7G"
    statement = ''' python %(project_src)s/UMI-tools/dedup_umi.py
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
@transform([dedup_alignments,indexMergedBAMs], 
           regex("(?:merged_)?(.+).bam(?:.bai)?"),
           r"\1.frag_length.tsv")
def getFragLengths(infile, outfile):
    ''' estimate fragment length distribution from read lengths'''

    intrack = re.match("(.+).bam(?:.bai)?", infile).groups()[0]

    statement = ''' python %(project_src)s/length_stats.py
                           -I %(intrack)s.bam
                           -S %(outfile)s
                           -L %(outfile)s.log
               '''
    if PARAMS["reads_paired"] == 1:
        statement += "--paired=%s" % (
            int(PARAMS["reads_length"]) - len(PARAMS["reads_bc_pattern"]))
    P.run()


###################################################################
@collate(getFragLengths,
         regex("(mapping|deduped).dir/.+\.frag_length.tsv"),
         r"\1.frag_lengths.load")
def loadFragLengths(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+\-.+\-.+).frag_length.tsv",
                         options=" -i Length")


###################################################################
@transform(dedup_alignments,
           suffix(".bam"), ".bam_stats.tsv")
def dedupedBamStats(infile, outfile):
    ''' Calculate statistics on the dedeupped bams '''

    statement = '''python %(scripts_dir)s/bam2stats.py
                         --force-output
                          < %(infile)s > %(outfile)s '''

    P.run()


###################################################################
@merge(dedupedBamStats, "deduped_bam_stats.load")
def loadDedupedBamStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).bam_stats.tsv")


###################################################################
@transform(dedup_alignments,
           suffix(".bam"),
           ".nspliced.txt")
def getNspliced(infile, outfile):
    ''' Calculate the number of reads spliced by grepping the
    cigar string for Ns'''

    statement = ''' samtools view %(infile)s
                  | awk '$6 ~ /N/'
                  | wc -l > %(outfile)s '''
    P.run()


###################################################################
@merge(getNspliced, "deduped_nspliced.load")
def loadNspliced(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).nspliced.txt",
                         cat="track",
                         has_titles=False,
                         header="track,nspliced",)


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
@follows(mkdir("saturation.dir"), run_mapping)
@subdivide(indexMergedBAMs, regex(".+/merged_(.+)\.[^\.]+\.bam.bai"),
       [r"saturation.dir/\1.%.3f.bam" % (1.0/(2 ** x))
        for x in range(1, 6)] +
       [r"saturation.dir/\1.%.3f.bam" % (x/10.0)
        for x in range(6, 11, 1)])
def subsetForSaturationAnalysis(infile, outfiles):
    '''Perform subsetting of the original BAM files, dedup and
    return the context stats. Test for resturn on investment for
    further sequencing of the same libraries '''

    track = re.match(".+/merged_(.+)\.[^\.]+\.bam", infile).groups()[0]
    infile = P.snip(infile, ".bai")
    statements = []
    statement_template = '''
                            python %%(project_src)s/UMI-tools/dedup_umi.py
                              %%(dedup_options)s
                              -I %(infile)s
                              -L %(outfile)s.log
                              -S %(outfile)s.tmp.bam
                              --subset=%(subset).3f;
                             checkpoint;
                             samtools sort %(outfile)s.tmp.bam %(outfile)s;
                             checkpoint;
                             samtools index %(outfile)s.bam '''
    for x in range(1, 6):
        subset = 1.0/(2 ** x)
        outfile = "saturation.dir/%%(track)s.%.3f" % subset
        statements.append(statement_template % locals())

    for x in range(6, 11):
        subset = x/10.0
        outfile = "saturation.dir/%%(track)s.%.3f" % subset
        statements.append(statement_template % locals())

    P.run()


###################################################################
@transform(subsetForSaturationAnalysis, suffix(".bam"), ".bamstats.tsv")
def subsetBamStats(infile, outfile):
    ''' Stats on the subset BAMs '''

    job_options = "-l mem_free=500M"
    statement = ''' python %(scripts_dir)s/bam2stats.py 
                    --force-output < %(infile)s > %(outfile)s '''
    P.run()


###################################################################
@merge(subsetBamStats, "subset_bam_stats.load")
def loadSubsetBamStats(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename= ".+/(.+-.+-.+)\.([0-9]+\.[0-9]+).bamstats.tsv",
                         cat="track,subset")


###################################################################
@transform([indexMergedBAMs, dedup_alignments],
           regex("(?:merged_)?(.+).bam(?:.bai)?"),
           add_inputs(generateContextBed),
           r"\1.reference_context.tsv")
def buildContextStats(infiles, outfile):
    ''' Find context from reads '''

    infile, reffile = infiles
    infile = re.match("(.+.bam)(?:.bai)?", infile).groups()[0]
    statement = ''' python %(scripts_dir)s/bam_vs_bed.py
                   --min-overlap=0.5
                   --log=%(outfile)s.log
                   %(infile)s %(reffile)s
                > %(outfile)s '''

    job_options = "-l mem_free=4G"
    P.run()


###################################################################
@collate(buildContextStats,
         regex("(mapping|deduped|saturation).dir/(?:[^/]+.dir/)?(.+).tsv"),
         r"\1_context_stats.load")
def loadContextStats(infiles, outfile):

    if "saturation" in infiles[0]:
        regex_filename = ".+/(.+-.+-.+)\.([0-9]+\.[0-9]+).reference_context.tsv"
        cat = "track,subset"
    else:
        regex_filename = ".+/(.+).reference_context.tsv"
        cat = "track"

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=regex_filename,
                         cat=cat)


###################################################################
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.flat.gtf.gz")
def flattenGeneSet(infile, outfile):
    ''' Get gtf with geneset flattened so that introns are definately
    intron sequence '''

    statement = '''python %(scriptsdir)s/gtf2gtf.py
                          --method=sort
                          --sort-order=gene+transcript
                          -I %(infile)s
                 | python %(scriptsdir)s/gtf2gtf.py
                          --method=merge-exons
                           -L %(outfile)s.log
                           -S %(outfile)s '''

    P.run()


###################################################################
@transform(dedup_alignments, suffix(".bam"),
           add_inputs(flattenGeneSet),
           ".splicing_index")
def calculateSplicingIndex(infiles, outfile):
    '''Calculate the splicing index: number of reads
    reads spliced at each junction divided by number of
    reads not spliced'''

    bamfile, gtffile = infiles
    PipelineiCLIP.calculateSplicingIndex(bamfile,
                                         gtffile,
                                         outfile,
                                         submit=True)


@merge(calculateSplicingIndex, "splicing_index.load")
def loadSplicingIndex(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="deduped.dir/(.+).splicing_index")

###################################################################
@follows(loadContextStats,
         
         loadDedupedBamStats,
         loadFragLengths,
         loadNspliced,
         loadDedupedUMIStats,
         loadSplicingIndex)
def MappingStats():
    pass

         
###################################################################
# Quality control and reproducibility
###################################################################
@follows(mkdir("reproducibility.dir"))
@collate(dedup_alignments, regex(".+/(.+\-.+)\-.+.bam"),
         r"reproducibility.dir/\1-agg.reproducibility.tsv.gz")
def calculateReproducibility(infiles, outfile):
    ''' Calculate cross-link reproducibility as defined by 
    Sugimoto et al, Genome Biology 2012 '''

    job_options = "-l mem_free=1G"
    infiles = " ".join(infiles)

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@follows(mkdir("reproducibility.dir"))
@merge(dedup_alignments,
       r"reproducibility.dir/agg-agg-agg.reproducibility.tsv.gz")
def reproducibilityAll(infiles, outfile):
    ''' Test weather sites from one file reproduce in other of different factors '''

    job_options = "-l mem_free=10G"
    infiles = " ".join(infiles)

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@follows(mkdir("reproducibility.dir"))
@transform(dedup_alignments,
           regex(".+/(.+).bam"),
           add_inputs("deduped.dir/%s*.bam" % PARAMS["experiment_input"]),
           r"reproducibility.dir/\1_vs_control.reproducibility.tsv.gz")
def reproducibilityVsControl(infiles, outfile):
    '''Test what fraction of the locations in each bam also appear
    in the control files'''

    track = infiles[0]
    if track in infiles[1:]:
        P.touch(outfile)
    else:
        job_options = "-l mem_free=1G"

        infiles = " ".join(infiles)

        statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                   -t %(track)s
                 | gzip > %(outfile)s '''
        P.run()


###################################################################
@merge(calculateReproducibility,
       r"reproducibility.dir/experiment_reproducibility.load")
def loadReproducibility(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile, cat="Experiment",
                         regex_filename=".+/(.+)-agg.reproducibility.tsv.gz",
                         options="-i Track -i fold -i level")


###################################################################
@transform(reproducibilityAll, regex("(.+)"),
           "reproducibility.dir/all_reproducibility.load")
def loadReproducibilityAll(infile, outfile):
    P.load(infile, outfile, "-i Track -i fold -i level")


###################################################################
@merge(reproducibilityVsControl,
       "reproducibility.dir/reproducibility_vs_control.load")
def loadReproducibilityVsControl(infiles, outfile):

    infiles = [infile for infile in infiles 
               if not PARAMS["experiment_input"] in infile]

    P.concatenateAndLoad(infiles, outfile, cat="Experiment",
                         regex_filename=".+/(.+\-.+)\-.+_vs_control.reproducibility.tsv.gz",
                         options = "-i File -i fold -i level")


###################################################################
@permutations(dedup_alignments, formatter(".+/(?P<TRACK>.+).bam"),
              2,
              "reproducibility.dir/{TRACK[0][0]}_vs_{TRACK[1][0]}.tsv.gz")
def computeDistances(infiles, outfile):
    ''' Compute the reproduciblity between each indevidual pair of samples
    this can then be readily converted to a distance measure'''

    track = infiles[0]
    infiles = " ".join(infiles)

    job_options="-l mem_free=2G"

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                   -t %(track)s
                   -m 1
                 | gzip > %(outfile)s '''

    P.run()


###################################################################
@merge(computeDistances,
       "reproducibility.dir/reproducibility_distance.load")
def loadDistances(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_vs_(.+).tsv.gz",
                         cat="File1,File2",
                         options="-i Track, -i File2")


###################################################################
@follows(loadReproducibility,
         loadReproducibilityAll,
         loadReproducibilityVsControl,
         loadDistances)
def reproducibility():
    pass


###################################################################
@follows(mkdir("counts.dir"))
@transform(dedup_alignments, regex(".+/(.+).bam"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"counts.dir/\1.tsv.gz")
def countReadsOverGenes(infiles, outfile):
    ''' use feature counts to quantify the number of tags on each gene'''

    bamfile, annotations = infiles

    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        nthreads=PARAMS['featurecounts_threads'],
        strand=1,
        options=PARAMS['featurecounts_options'] + ' -O')


###################################################################
@merge(countReadsOverGenes,
       "counts.dir/track_counts.tsv.gz")
def mergeCounts(infiles, outfile):
    '''Merge feature counts data into one table'''

    infiles = " ".join(infiles)

    statement=''' python %(scriptsdir)s/combine_tables.py
                         -c 1
                         -k 7
                         --regex-filename='(.+).tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(mergeCounts, suffix(".tsv.gz"), ".load")
def loadCounts(infile, outfile):

    P.load(infile, outfile, options="-i geneid")

###################################################################
# Analysis
###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform(dedup_alignments, regex(".+/(.+).bam"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"gene_profiles.dir/\1.tsv")
def calculateGeneProfiles(infiles, outfile):
    ''' Calculate metagene profiles over protein coding genes
    for each sample'''

    infile, reffile = infiles
    statement = '''python %(scripts_dir)s/bam2geneprofile.py
                           --method=geneprofilewithintrons
                           --bam-file=%(infile)s
                           --gtf-file=%(reffile)s
                           --normalize-transcript=total-sum
                           --normalize-profile=area
                           --log=%(outfile)s.log
                           --output-filename-pattern=%(outfile)s.%%s
                           > %(outfile)s '''

    P.run()


@merge(calculateGeneProfiles,
       "gene_profiles.load")
def loadGeneProfiles(infiles, outfile):

    infiles = [infile + ".geneprofilewithintrons.matrix.tsv.gz"
               for infile in infiles]

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename='.+/(.+)\-(.+)\-(.+).tsv.geneprofilewithintrons.matrix.tsv.gz',
                          cat = "factor,condition,rep",
                          options = "-i factor -i condition -i rep")

###################################################################
@follows(mapping_qc)
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.exons.gtf.gz")
def transcripts2Exons(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons
    only keeps those exons longer than 100 base pairs'''
    
    tmp_outfile = P.snip(outfile, ".gtf.gz") + ".tmp.gtf.gz"
    PipelineiCLIP.removeFirstAndLastExon(infile, tmp_outfile)

    statement = '''python %(scripts_dir)s/gff2bed.py 
                          --is-gtf 
                           -I %(tmp_outfile)s 
                           -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
                  | mergeBed -i stdin -s -d 100 -c 4 -o distinct
                  | awk '($3-$2) > 100 {print}'
                  | awk 'BEGIN{OFS="\\t"} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py --as-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()

    os.unlink(tmp_outfile)


###################################################################
@follows(mapping_qc)
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.introns.gtf.gz")
def transcripts2Introns(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons'''

    tmp_outfile = P.snip(outfile, ".gtf.gz") + ".tmp.gtf.gz"
    PipelineiCLIP.removeFirstAndLastExon(infile, tmp_outfile)

    statement = '''python %(scripts_dir)s/gtf2gtf.py
                           -I %(tmp_outfile)s
                          --log=%(outfile)s.log
                           --method=exons2introns
                  | python %(scripts_dir)s/gff2bed.py
                          --is-gtf
                           -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
                  | mergeBed -i stdin -s -d 100 -c 4 -o distinct
                  | awk 'BEGIN{OFS="\\t"} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py
                          --as-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()

    os.unlink(tmp_outfile)

###################################################################
@product(dedup_alignments,
         formatter(".+/(?P<TRACK>.+).bam"),
         [transcripts2Exons, transcripts2Introns],
         formatter(".+/geneset_all.(?P<INTERVALTYPE>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{INTERVALTYPE[1][0]}.log")
def calculateExonProfiles(infiles, outfile):

    infile, reffile = infiles

    outfile = P.snip(outfile, ".log")
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=intervalprofile
                          --bam-file=%(infile)s
                          --gtf-file=%(reffile)s
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --resolution-upstream=50
                          --resolution-downstream=50
                          --extension-upstream=50
                          --extension-downstream=50
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.log '''

    P.run()

@merge(calculateExonProfiles,
       "exon_profiles.load")
def loadExonProfiles(infiles, outfile):

    infiles = [P.snip(infile,".log") + ".intervalprofile.matrix.tsv.gz"
               for infile in infiles]

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename='.+/(.+)\-(.+)\-(.+).(exons|introns).intervalprofile.matrix.tsv.gz',
                          cat = "factor,condition,rep,interval",
                          options = "-i factor -i condition -i rep")
###################################################################
@product(dedup_alignments,
         formatter(".+/(?P<TRACK>.+).bam"),
         [transcripts2Exons, transcripts2Introns],
         formatter(".+/refcoding.(?P<INTERVALTYPE>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{INTERVALTYPE[1][0]}.tssprofile.log")
def calculateExonTSSProfiles(infiles, outfile):

    infile, reffile = infiles

    outfile = P.snip(outfile, ".tssprofile.log")
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=tssprofile
                          --bam-file=%(infile)s
                          --gtf-file=%(reffile)s
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --resolution-upstream=100
                          --resolution-downstream=100
                          --extension-inward=100
                          --extension-outward=100
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.tssprofile.log '''

    P.run()


###################################################################
@follows(calculateExonTSSProfiles,
         calculateGeneProfiles,
         calculateExonProfiles,
         loadExonProfiles,
         loadGeneProfiles)
def profiles():
    pass


###################################################################
# Calling significant clusters
###################################################################
@follows(mkdir("clusters.dir"), mapping_qc)
@subdivide(dedup_alignments,
           regex(".+/(.+).bam"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           [r"clusters.dir/\1.bg.gz",
            r"clusters.dir/\1.bed.gz"])
def callSignificantClusters(infiles, outfiles):
    '''Call bases as significant based on mapping depth in window
    around base'''

    bam, gtf = infiles
    PipelineiCLIP.callClusters(bam, gtf, outfiles)


###################################################################
@collate(callSignificantClusters,
         regex("(.+/.+\-.+)\-(.+)\.bed.gz"),
         r"\1.reproducible.bed.gz")
def callReproducibleClusters(infiles, outfile):
    '''Find clusters that appear in more than one replicate'''

    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           PARAMS["clusters_min_reproducible"])


###################################################################
@transform(callSignificantClusters, suffix(".bg.gz"), ".count_bases")
def countCrosslinkedBases(infile, outfile):
    ''' Count number of crosslinked bases within gene models'''

    statement = ''' zcat %(infile)s |
                    wc -l > %(outfile)s '''

    P.run()


###################################################################
@merge(countCrosslinkedBases, "cross_linked_bases.load")
def loadCrosslinkedBasesCount(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).count_bases",
                         header="track,count",
                         cat="track",
                         has_titles=False)


###################################################################
@transform([callSignificantClusters, callReproducibleClusters],
           suffix(".bed.gz"),
           r".cluster_count")
def countClusters(infile, outfile):
    '''Count the number of significant clusters'''

    statement = '''zcat %(infile)s | wc -l > %(outfile)s'''
    P.run()


###################################################################
@merge(countClusters, "cluster_counts.load")
def loadClusterCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).(R[0-9]+|reproducible).cluster_count",
                         header="sample,replicate,count",
                         cat="sample,replicate",
                         has_titles=False)


###################################################################
@transform([callSignificantClusters, callReproducibleClusters],
           suffix(".bed.gz"),
           add_inputs(generateContextBed),
           ".context_stats.tsv.gz")
def getClusterContextStats(infiles, outfile):
    '''Generate context stats for called clusters'''

    clusters, context = infiles
    tmp = P.getTempFilename()

    statement = '''  zcat %(clusters)s | sort -k1,1 -k2,2n > %(tmp)s.bed;
                     checkpoint;
                     python %(scriptsdir)s/bam_vs_bed.py
                           -a %(tmp)s.bed
                           -b %(context)s -S  %(outfile)s;
                     checkpoint;
                     rm %(tmp)s.bed'''

    P.run()


###################################################################
@merge(getClusterContextStats, "clusters.dir/cluster_context_stats.load")
def loadClusterContextStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=
                         "clusters.dir/(.+).context_stats.tsv.gz")


###################################################################
@follows(callSignificantClusters,
         loadCrosslinkedBasesCount,
         loadClusterCounts,
         loadClusterContextStats)
def clusters():
    pass


###################################################################
# Motifs
###################################################################

@transform([callSignificantClusters, callReproducibleClusters],
           regex(".+/(.+).bed.gz"),
           r"clusters.dir/\1.fa")
def clusters2fasta(infile, outfile):
    '''convert clusters to fasta ready for motif calling
       only keep those greater than 8 bp long'''

    statement = '''
               zcat %(infile)s 
              | awk '$3-$2 > 8'   
              | python %(scriptsdir)s/bed2fasta.py
                       -g %(genome_dir)s/%(genome)s
                       -m dustmasker
                       --use-strand
                       --output-mode=segments
                       
                       -L %(outfile)s.log
              | sed 's/[ |\:]/_/g' > %(outfile)s
              '''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.fa.gz")
def getReferenceGenesetFasta(infile, outfile):
    '''Collapse genesets onto single intervals and get fasta'''

    logfile = P.snip(outfile, ".fa.gz") + ".log"

    statement = '''python %(scriptsdir)s/gtf2gtf.py
                           --method=merge-transcripts
                            -I %(infile)s -L %(logfile)s
                 | python %(scriptsdir)s/gff2bed.py 
                            -L %(logfile)s
                 | python %(scriptsdir)s/bed2fasta.py
                            -g %(genome_dir)s/%(genome)s
                            -m dustmasker
                            --use-strand
                            -L %(outfile)s.log 
                 | gzip > %(outfile)s'''

    P.run()


###################################################################
@transform(getReferenceGenesetFasta,
           regex(".+/(.+).fa.gz"),
           r"meme.dir/\1.model")
def getMEMEBackgroundModel(infile, outfile):
    '''Get a background markov model for MEME based on the input
    sequence to the cluster calling algorithmn '''

    statement = '''zcat %(infile)s
                 | fasta-get-markov -m %(meme_background_order)s
                                    --norc
                 > %(outfile)s 2> %(outfile)s.log '''
    P.run()


###################################################################
@follows(mkdir("meme.dir"))
@transform(clusters2fasta, regex(".+/(.+).fa"),
           r"meme.dir/\1.meme")
def runMeme(infiles, outfile):
    '''Run Meme to find motifs. All intervals are currently used
    all of each interval is used'''

    foreground = infiles
    PARAMS["meme_revcomp"] = False
    PipelineMotifs.PARAMS = PARAMS
    
    tmpfile = P.getTempFilename(shared=True)
    logfile = outfile + ".log"
    PipelineiCLIP.subsampleNReadsFromFasta(foreground, tmpfile,
                                           PARAMS["meme_max_sequences"],
                                           logfile)
    PipelineMotifs.runMEMEOnSequences(tmpfile, outfile)

    os.unlink(tmpfile)


 ###################################################################
@merge(runMeme, "meme_summary.load")
def loadMemeSummary(infiles, outfile):
    ''' load information about motifs into database '''

    outf = P.getTempFile(".")
    outf.write("track\n")
    
    for infile in infiles:
        if IOTools.isEmpty(infile):
            continue
        motif = P.snip(infile, ".meme")
        outf.write("%s\n" % motif)

    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)


###################################################################
@follows(mkdir("dreme.dir"))
@transform(clusters2fasta, regex(".+/(.+).fa"),
           r"dreme.dir/\1.txt")
def runDREME(infile, outfile):
    '''run Dreme on full set of clusters, using shuffled
    input as background'''

    PipelineMotifs.runDREME(infile, outfile)

###################################################################
@follows(loadMemeSummary)
def meme():
    pass


###################################################################
@follows(meme, runDREME)
def motifs():
    pass



###################################################################
# Data Export
###################################################################
@collate(dedup_alignments,
         regex("(.+\-.+)\-(.+).bam"),
         r"\1.union.bam")
def makeUnionBams(infiles, outfile):
    '''Merge replicates together'''

    outfile = os.path.abspath(outfile)

    if len(infiles) == 1:
        infile = os.path.abspath(infiles[0])
        statement = '''ln -sf %(infile)s %(outfile)s;
                       checkout;
 
                       ln -sf %(infile)s.bai %(outfile)s.bai;'''
    else:

        statement = ''' samtools merge -f %(outfile)s %(infiles)s;
                        checkpoint;

                        samtools index %(outfile)s'''

    infiles = " ".join(infiles)

    P.run()

@follows(mkdir("bigWig"))
@subdivide([dedup_alignments,
            makeUnionBams], 
           regex(".+/(.+).bam"),
          [r"bigWig/\1_plus.bw",
           r"bigWig/\1_minus.bw"])
def generateBigWigs(infile, outfiles):
    '''Generate plus and minus strand bigWigs from BAM files '''

    out_pattern = P.snip(outfiles[0], "_plus.bw")
    statement = '''python %(project_src)s/iCLIP2bigWig.py
                          -I %(infile)s
                          -L %(out_pattern)s.log
                          %(out_pattern)s '''

    P.run()

###################################################################
@follows(mkdir("export/hg19"))
@transform(generateBigWigs,
           regex("bigWig/(.+)"),
           r"export/hg19/\1")
def linkBigWig(infile, outfile):
    '''Link bigwig files to export directory'''
    
    try:
        os.symlink(os.path.abspath(infile), os.path.abspath(outfile))
    except OSError:
        os.unlink(outfile)
        os.symlink(os.path.abspath(infile), os.path.abspath(outfile))


###################################################################
@merge(linkBigWig, "export/hg19/tagwig_trackDb.txt")
def generateBigWigUCSCFile(infiles, outfile):
    '''Generate track configuration for exporting wig files '''


    track_template = '''
          track tagwig_%(track)s_%(strand)s
          parent tagwig_%(track)s
          bigDataUrl %(track_data_URL)s
          shortLabel %(short_label)s
          longLabel %(long_label)s
          color %(color)s
          type bigWig %(negate)s'''

    overlap_template = '''
       track tagwig_%(track)s
       parent clipwig
       shortLabel %(short_label)s
       longLabel %(long_label)s
       autoScale on
       visibility full
       container multiWig
       type bigWig
       aggregate solidOverlay
       maxHeightPixels 16:16:32
       alwaysZero on'''

    stanzas = {}
    for infile in infiles:
        track, strand = re.match(
            ".+/(.+-.+)_(plus|minus).bw", infile).groups()

        negate = ""

        if strand == "minus":
            negate = '''
          negateValues on'''
            color="255,0,0"
        else:
            color="0,0,255"
            
        track_data_URL = os.path.basename(infile)
        short_label = track + "iCLIP tags"
        long_label = "iCLIP tags from track %s" \
                     % track
        
        if track not in stanzas:
            stanzas[track] = overlap_template % locals()

        stanzas[track] += "\n" + track_template % locals()

    composite_stanaz = '''
    track clipwig
    shortLabel iCLIP tags
    longLabel Raw iCLIP tags
    superTrack on
    alwaysZero on
    maxHeightPixels 16:16:32'''

    output = "\n".join([composite_stanaz] + list(stanzas.values()))

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(output)


###################################################################
@follows(mkdir("export/hg19"))
@transform([callSignificantClusters, callReproducibleClusters],
           regex("clusters.dir/(.+).bed.gz"),
           r"export/hg19/\1.bigBed")
def exportClusters(infile, outfile):
    ''' Add a track line to cluster files and export to export dir '''
   
    PipelineiCLIP.clustersToBigBed(infile, outfile)


###################################################################
@merge(exportClusters, "export/hg19/clusters_trackDb.txt")
def generateClustersUCSC(infiles, outfile):

    PipelineiCLIP.makeClustersUCSC(infiles, outfile, "pipelineClusters",
                                   "Clusters from iCLIP pipeline")


###################################################################
@merge([generateClustersUCSC, generateBigWigUCSCFile],
       "export/hg19/trackDb.txt")
def mergeTrackDbs(infiles, outfile):

    to_cluster = False
    infiles = " ".join(infiles)
    statement = "cat %(infiles)s > %(outfile)s"
    P.run()


###################################################################
@follows(mkdir("export"))
@originate(["export/hub.txt",
            "export/genomes.txt"])
def makeHubFiles(outfiles):

    hub_file = '''
    hub iCLIPPipeline%(version)s
    shortLabel CGAT iCLIP Pipelines
    longLabel All browser tracks CGAT iCLIP pipeline run
    genomesFile genomes.txt
    email i.sudbery@sheffield.ac.uk''' % PARAMS

    with IOTools.openFile("export/hub.txt", "w") as outf:
        outf.write(hub_file)

    genomes_file = '''
    genome hg19
    trackDb hg19/trackDb.txt'''

    with IOTools.openFile("export/genomes.txt", "w") as outf:
        outf.write(genomes_file)


###################################################################
@follows(mergeTrackDbs, makeHubFiles)
def export():
    pass


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, mapping, MappingStats, reproducibility,
         profiles, clusters, motifs, export)
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
    P.run_report( clean = True )


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
