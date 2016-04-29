##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
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
"""===========================
Pipeline template
===========================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_umi_paper.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
import sys
import os
import urllib
import sqlite3
from bs4 import BeautifulSoup
import pandas as pd
from rpy2.robjects import r as R

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.Sra as SRA
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator

import CGATPipelines.PipelinePreprocess as PipelinePreprocess
import CGATPipelines.PipelineMapping as PipelineMapping
import PipelineScRNASeq

###################################################
# Pipeline configuration
###################################################
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False},
    only_import=__name__ != "__main__")

PARAMS = P.PARAMS


def makeSoup(address):
    sock = urllib.urlopen(address)
    htmlSource = sock.read()
    soup = BeautifulSoup(htmlSource)
    return soup


###############################################################################
# Section - START - GSE53638 - Soumillon et al 2014
###############################################################################

@follows(mkdir("GSE53638"))
@originate(["GSE53638/SRR1058003.sra",
            "GSE53638/SRR1058023.sra",
            "GSE53638/SRR1058032.sra",
            "GSE53638/SRR1058038.sra"])
def downloadGGSE53638(outfile):
    ''' download the sra files '''

    address_base = 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP034%2FSRP034712'

    sra = os.path.basename(outfile).replace(".sra", "")

    statement = '''cd GSE53638;
    wget %(address_base)s/%(sra)s/%(sra)s.sra'''

    P.run()


@mkdir("GSE53638/fastqs.dir")
@transform(downloadGGSE53638,
           regex("GSE53638/(\S+).sra"),
           r"GSE53638/fastqs.dir/\1_1.fastq.gz")
def extractGGSE53638(infile, outfile):
    ''' extract the fastqs from the SRA '''

    statement = SRA.extract(infile, "GSE53638/fastqs.dir")

    P.run()


@subdivide(extractGGSE53638,
           regex("(\S+)_1.fastq.gz"),
           r"\1_UMI_*.fastq.gz")
def extractUMIsAndFilterGSE53638(infile, outfiles):
    ''' extract UMIs from read 1 and filter as per Soumillon et al 2014'''

    UMI_fastq = infile
    fastq = infile.replace("_1.fastq.gz", "_2.fastq.gz")

    # ToDo: replace with PARAMS['soumillon_barcodes']
    barcode_inf = "Soulmillon_et_al_2014_barcodes.tsv"
    barcode_inf = IOTools.openFile(barcode_inf, "r")
    barcodes = []

    # different barcodes for Differentiation 1 and Differentiation 3
    if "SRR1058003" in UMI_fastq or "SRR1058023" in UMI_fastq:
        start = -1

    elif "SRR1058032" in UMI_fastq or "SRR1058038" in UMI_fastq:
        start = -2

    for line_number, line in enumerate(barcode_inf.read().splitlines(), start):
        if line_number % 3 == 0:
            barcodes.append(line)

    PipelineScRNASeq.extractUMIsAndFilterFastq(fastq, UMI_fastq, barcodes,
                                               submit=True)


@originate("GSE53638/ERCC.fasta")
def createERCCFastaGSE53638(outfile):
    ''' Create fasta file from ERCC spike in txt file'''

    outf = IOTools.openFile(outfile, "w")

    # ToDo - where the infile coming from? Parameterise infile location
    with IOTools.openFile("data/cms_095047.txt", "r") as inf:
        next(inf)
        for line in inf:
            id = ">" + " ". join(line.split("\t")[:-1])
            seq = line.split("\t")[-1]
            split_at = 60
            seq = [seq[i:i+split_at] for i in range(0, len(seq), split_at)]
            outf.write("%s\n%s" % (id, "\n".join(seq)))


@follows(mkdir("GSE53638/geneset.dir"))
@merge("/ifs/mirror/annotations/hg19_ensembl75/geneset_all.gtf.gz",
       "GSE53638/geneset.dir/reference.gtf.gz")
def buildReferenceGeneSetGSE53638(infile, outfile):
    '''filter to protein-coding transcripts and merge transcripts
    '''

    # Returns an sqlite3 database handle.
    # ToDo: PARAMS['soumillon_annotaions_database']
    dbh = sqlite3.connect("/ifs/mirror/annotations/hg19_ensembl75/csvdb")

    select_cmd = """
    SELECT DISTINCT transcript_id from transcript_info
    WHERE gene_biotype = 'protein_coding'
    AND transcript_biotype = 'protein_coding'
    """ % locals()

    select = dbh.execute(select_cmd)

    tmp_out = P.getTempFilename(shared=True)

    with IOTools.openFile(tmp_out, "w") as outf:
        outf.write("\n".join((x[0] for x in select)) + "\n")
    outf.close()

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py
    --method=filter
    --filter-method=transcript
    --map-tsv-file=%(tmp_out)s
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
    --method=merge-exons
    --merge-exons-distance=10
    --with-utr
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
    --method=set-transcript-to-gene
    | gzip
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmp_out)


@transform(buildReferenceGeneSetGSE53638,
           suffix(".gtf.gz"),
           ".fa")
def buildReferenceTranscriptomeGSE53638(infile, outfile):
    '''build reference transcriptome'''

    gtf_file = P.snip(infile, ".gz")

    # ToDo: os.path.join(
    # PARAMS['soumillon_genome_directory'],
    # PARAMS['soumillon_genome'] + ".fasta")
    genome_file = "/ifs/mirror/genomes/plain/hg19.fasta"

    statement = '''
    zcat %(infile)s
    | awk '$3~/exon|UTR/' > %(gtf_file)s;
    gtf_to_fasta %(gtf_file)s %(genome_file)s %(outfile)s;
    checkpoint;
    samtools faidx %(outfile)s
    '''
    P.run()

    dest = P.snip(os.path.abspath(gtf_file), ".gtf") + ".gff"
    if not os.path.exists(dest):
        os.symlink(os.path.abspath(gtf_file), dest)


@merge([createERCCFastaGSE53638,
        buildReferenceTranscriptomeGSE53638],
       "GSE53638/geneset.dir/refcoding_plus_ERCC.fa")
def mergeReferenceTranscriptomeERCCGSE53638(infiles, outfile):
    ''' merge fasta from reference transcriptome with ERCC spike in fasta and mitochondrial DNA'''

    mtGenome = P.getTempFilename()

    # ToDo: os.path.join(
    # PARAMS['soumillon_genome_directory'],
    # PARAMS['soumillon_genome'] + ".fasta")
    
    # Extract Mt chromosome
    FastaIt = FastaIterator.FastaIterator(
        IOTools.openFile("/ifs/mirror/genomes/plain/hg19.fasta", "r"),
        split=60)

    with IOTools.openFile(mtGenome, "w") as outf:
        for entry in FastaIt:
            if entry.title == "chrM":
                outf.write("%s\n" % entry)

    ERCC, reference = infiles

    # sed command ensures fasta entries are correctly labelled
    statement = '''
    cat %(ERCC)s %(reference)s %(mtGenome)s |
    sed 's/>[0-9]\+ />/g' > %(outfile)s;
    samtools faidx %(outfile)s'''
    P.run()


@transform(mergeReferenceTranscriptomeERCCGSE53638,
           suffix(".fa"),
           ".sa")
def indexMergedFastaGSE53638(infile, outfile):
    ''' build BWA index for merged fasta'''

    prefix = P.snip(outfile, ".sa")

    # build raw index
    statement = '''
    bwa index %(infile)s -p %(prefix)s >> %(outfile)s.log 2>&1
    '''
    P.run()


@mkdir("GSE53638/transcriptome.dir")
@follows(indexMergedFastaGSE53638)
@transform(extractUMIsAndFilterGSE53638,
           regex("GSE53638/fastqs.dir/(\S+)_UMI_(\S+).fastq.gz"),
           add_inputs(indexMergedFastaGSE53638),
           r"GSE53638/transcriptome.dir/\1_UMI_\2.trans.bam")
def mapBWAAgainstGenesetGSE53638(infiles, outfile):
    ''' map reads using BWA against transcriptome data

    bwa parameterised according to soumillon et al 2014:
    -l 24 = seed length - 24 bp
    -k 2 = default number of mismatches allowed in seed - 2
    -n 0.04 = default percentage of mismatches allowed across read - 4%

    non-unique alignments will NOT be removed from the final bam
    '''

    infile, reference = infiles
    job_threads = 2
    job_options = "-l mem_free=1.9G"
    bwa_aln_options = "-l 24 -k 2 -n 0.04"
    bwa_index_dir = os.path.dirname(reference)
    genome = P.snip(os.path.basename(reference), ".sa")
    bwa_threads = job_threads
    bwa_samse_options = ""
    m = PipelineMapping.BWA(remove_non_unique=0,
                            strip_sequence=0,
                            set_nh=1)

    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("GSE53638/dedup_unique.dir"),
         mkdir("GSE53638/dedup_percentile.dir"),
         mkdir("GSE53638/dedup_cluster.dir"),
         mkdir("GSE53638/dedup_adjacency.dir"),
         mkdir("GSE53638/dedup_directional_adjacency.dir"),
         mapBWAAgainstGenesetGSE53638)
@subdivide(mapBWAAgainstGenesetGSE53638,
           regex("GSE53638/transcriptome.dir/(\S+)_UMI_(\S+).trans.bam"),
           [r"GSE53638/dedup_unique.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE53638/dedup_percentile.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE53638/dedup_cluster.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE53638/dedup_adjacency.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE53638/dedup_directional_adjacency.dir/\1_UMI_\2_deduped.trans.bam"])
def dedupGSE53638(infile, outfiles):
    ''' perform deduping with various methods'''

    for outfile in outfiles:
        outfile_tmp = outfile.replace("_deduped.trans.bam", "_temp.trans.bam")
        outfile_snip = P.snip(outfile, ".bam")

        method = P.snip(os.path.basename(
            os.path.dirname(outfile).replace("dedup_", "")), ".dir")
        
        if method == "directional_adjacency":
            method = "directional-adjacency"
        
        if method == "cluster":
            options = "--further-stats"
        else:
            options = ""

        statement = '''
        umi_tools dedup
        --per-contig --method=%(method)s
        -I %(infile)s -S %(outfile_tmp)s -L %(outfile)s.log
        --output-stats=%(outfile)s.stats %(options)s;
        checkpoint;
        samtools sort %(outfile_tmp)s %(outfile_snip)s;
        checkpoint;
        samtools index %(outfile)s;
        checkpoint;
        rm -rf %(outfile_tmp)s'''
        P.run()


@collate([mapBWAAgainstGenesetGSE53638,
          dedupGSE53638],
         regex("GSE53638/(\S+).dir/(\S+)_UMI_(\S+).bam"),
         r"GSE53638/\1.dir/\2_UMI_\3.gene.counts.tsv")
def countGenesGSE53638(infile, outfile):
    ''' summarise counts per gene '''

    job_memory = "2G"
    PipelineScRNASeq.countAlignmentsPerGene(
        infile[0], outfile, mapq_threshold=10, submit=True)


@follows(mkdir("GSE53638/figures.dir"))
@collate(countGenesGSE53638,
         regex("GSE53638/(\S+).dir/(\S+)_UMI_(\S+).trans.gene.counts.tsv"),
         r"GSE53638/figures.dir/\1_\2_merged_gene_counts.tsv")
def mergeCountsGSE53638(infiles, outfile):
    ''' merge gene counts for each cell in each sample into a single
    table per sample'''

    job_memory = "4G"
    PipelineScRNASeq.getGeneCounts(infiles, outfile, submit=True)


@collate(dedupGSE53638,
         regex(r"GSE53638/(\S+)/(\S+)_UMI_\S+_deduped.trans.bam"),
         r"GSE53638/\1/\2_edit_distances.tsv")
def editDistanceDedupGSE53638(infiles, outfile):
    '''concatenate edit distance stats'''

    infiles = [infile + ".stats_edit_distance.tsv" for infile in infiles]

    job_memory = "1G"
    PipelineScRNASeq.summariseEditDistances(infiles, outfile, submit=True)


@follows(mkdir("GSE53638/figures.dir"))
@collate(editDistanceDedupGSE53638,
         regex("GSE53638/(\S+)/(\S+)_edit_distances.tsv"),
         [r"GSE53638/figures.dir/\2_merged_edit_distances.tsv",
          r"GSE53638/figures.dir/\2_edit_distances.png"])
def mergeAndPlotEditDistancesGSE53638(infiles, outfiles):
    '''merge edit distance stats for all deduplication methods and plot'''

    outfile, plot_out = outfiles

    job_memory = "1G"
    PipelineScRNASeq.mergeAndPlotEditDistances(
        infiles, outfile, plot_out, submit=True)


@collate(mergeCountsGSE53638,
         regex("GSE53638/figures.dir/dedup_(\S+)_(\S+)_merged_gene_counts.tsv"),
         [r"GSE53638/figures.dir/\2_ERCC_sd.png",
          r"GSE53638/figures.dir/\2_ERCC_mean.png"])
def plotERCCAccuracyGSE53638(infiles, outfiles):
    ''' merge the counts tables per sample, select the ERCC spike ins and
    plot'''

    PipelineScRNASeq.plotERCCAccuracy(infiles, outfiles, submit=True)


@merge(mergeCountsGSE53638,
       "GSE53638/figures.dir/heatmap.log")
def plotHeatmapsGSE53638(infiles, outfile):
    ''' plot heatmaps for each dedup method'''

    PipelineScRNASeq.makeHeatmapsAndPCA(infiles, outfile, submit=True)


@transform(plotHeatmapsGSE53638,
           regex("GSE53638/figures.dir/heatmap.log"),
           r"GSE53638/figures.dir/variance.png")
def plotVarianceGSE53638(infile, outfile):
    ''' plot Variance explained in first x PCs each dedup method'''

    PipelineScRNASeq.plotVarianceGSE53638(
        infile, outfile, submit=True)


@transform(plotHeatmapsGSE53638,
           regex("GSE53638/figures.dir/heatmap.log"),
           "GSE53638/figures.dir/vector_expression.png")
def plotVectorsVsExpressionGSE53638(infile, outfile):
    ''' plot Vectors vs Expression for first x PCs for each dedup method '''

    PipelineScRNASeq.plotVectorsVsExpressionGSE53638(
        infile, outfile, submit=True)


@collate(mergeCountsGSE53638,
         regex("GSE53638/figures.dir/(\S+)_SRR(\S+)_merged_gene_counts.tsv"),
         r"GSE53638/figures.dir/SRR\2_cv.png")
def plotCVGSE53638(infiles, plotfile):
    ''' Calculate CV for each method and plot'''

    job_memory = "4G"
    normalise_method = "total-count"

    PipelineScRNASeq.plotCV(
        infiles, plotfile, normalise_method,
        submit=True, job_memory=job_memory)


@follows(countGenesGSE53638,
         mergeAndPlotEditDistancesGSE53638,
         mergeCountsGSE53638,
         plotERCCAccuracyGSE53638,
         plotHeatmapsGSE53638,
         plotCVGSE53638,
         plotVarianceGSE53638,
         plotVectorsVsExpressionGSE53638)
def GSE53638():
    pass

###############################################################################
# Section END
###############################################################################

###############################################################################
# Section - START - GSE65525 - Klein et al 2015
###############################################################################


@follows(mkdir("GSE65525"))
@originate(["GSE65525/SRR1784310.sra",
            "GSE65525/SRR1784313.sra",
            "GSE65525/SRR1784314.sra",
            "GSE65525/SRR1784315.sra"])
def downloadGSE65525(outfile):
    ''' download SRAs '''

    base = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP053%2FSRP053052"

    sra = P.snip(os.path.basename(outfile), ".sra")

    statement = '''cd GSE65525; wget %(base_address)s/%(sra)s/%(sra)s.sra'''

    P.run()


@mkdir("GSE65525/fastqs.dir")
@transform(downloadGSE65525,
           regex("GSE65525/(\S+).sra"),
           r"GSE65525/fastqs.dir/\1_1.fastq.gz")
def extractGSE65525(infile, outfile):
    ''' extract fastqs '''
    statement = SRA.extract(infile, "GSE65525/fastqs.dir")

    P.run()


@subdivide(extractGSE65525,
           regex("GSE65525/fastqs.dir/(\S+)_1.fastq.gz"),
           r"GSE65525/fastqs.dir/\1_UMI_*.fastq.gz")
def extractUMIsAndFilterGSE65525(infile, outfiles):
    ''' extract UMIs from read 1 and filter as per Allon et al 2015'''

    # UMIs are in read1, sequence for alignment in read 2
    UMI_fastq = infile
    fastq = infile.replace("_1.fastq.gz", "_2.fastq.gz")

    barcodes1_infile = "data/gel_barcode1_list.txt"
    barcodes2_infile = "data/gel_barcode2_list.txt"

    sample = P.snip(os.path.basename(infile), "_1.fastq.gz")

    # These are the number of cell barcodes according to Klein et al
    sample2cellbarcodes = {
        "SRR1784310": 935,
        "SRR1784313": 301,
        "SRR1784314": 682,
        "SRR1784315": 799
    }

    cell_barcodes = sample2cellbarcodes[sample]

    job_memory = "2G"

    PipelineScRNASeq.extractUMIsAndFilterFastqGSE65525(
        UMI_fastq, fastq, barcodes1_infile, barcodes2_infile, cell_barcodes,
        submit=True)


@mkdir("GSE65525/processed.dir")
@transform(extractUMIsAndFilterGSE65525,
           regex("GSE65525/fastqs.dir/(\S+)_UMI_(\S+).fastq.gz"),
           r"GSE65525/processed.dir/trimmed-\1_UMI_\2.fastq.gz")
def processReadsGSE65525(infile, outfile):
    ''' process the reads with trimmomatic as per Klein et al 2015 '''

    track = P.snip(os.path.basename(infile), ".fastq.gz")

    threads = 1
    job_memory = "7G"

    # as per Allon et al 2015
    trimmomatic_options = "LEADING:28 SLIDINGWINDOW:4:20 MINLEN:19"

    m = PipelinePreprocess.MasterProcessor(
        threads=threads)

    m.add(PipelinePreprocess.Trimmomatic(
        trimmomatic_options, threads=threads))

    statement = m.build((infile,), "GSE65525/processed.dir/trimmed-", track)

    P.run()


@follows(mkdir("GSE65525/geneset.dir"))
@merge("/ifs/mirror/annotations/mm10_ensembl78/ensembl.dir/geneset_all.gtf.gz",
       "GSE65525/geneset.dir/reference.gtf.gz")
def buildReferenceGeneSetGSE65525(infile, outfile):
    '''filter to protein-coding transcripts and merge transcripts
    '''

    # Returns an sqlite3 database handle
    # ToDo PARAMS
    dbh = sqlite3.connect(
        "/ifs/mirror/annotations/mm10_ensembl78/csvdb")

    select_cmd = """
    SELECT DISTINCT transcript_id from transcript_info
    WHERE gene_biotype = 'protein_coding'
    AND transcript_biotype = 'protein_coding'
    """ % locals()

    select = dbh.execute(select_cmd)

    tmp_out = P.getTempFilename(shared=True)

    with IOTools.openFile(tmp_out, "w") as outf:
        outf.write("\n".join((x[0] for x in select)) + "\n")
    outf.close()

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py
    --method=filter
    --filter-method=transcript
    --map-tsv-file=%(tmp_out)s
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
    --method=merge-exons
    --merge-exons-distance=10
    --with-utr
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
    --method=set-transcript-to-gene
    | gzip
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmp_out)


@transform(buildReferenceGeneSetGSE65525,
           suffix(".gtf.gz"),
           ".fa")
def buildReferenceTranscriptomeGSE65525(infile, outfile):
    '''build reference transcriptome'''

    gtf_file = P.snip(infile, ".gz")

    # ToDo # PARAMS
    genome_file = "/ifs/mirror/genomes/plain/mm10.fasta"

    # sed command ensures fasta entries are correctly labelled
    statement = '''
    zcat %(infile)s
    | awk '$3~/exon|UTR/' > %(gtf_file)s;
    gtf_to_fasta %(gtf_file)s %(genome_file)s %(outfile)s;
    checkpoint;
    sed -i 's/>[0-9]\+ />/g' %(outfile)s;
    checkpoint;
    samtools faidx %(outfile)s
    '''
    P.run()

    dest = P.snip(os.path.abspath(gtf_file), ".gtf") + ".gff"
    if not os.path.exists(dest):
        os.symlink(os.path.abspath(gtf_file), dest)


@transform(buildReferenceTranscriptomeGSE65525,
           suffix(".fa"),
           ".1.ebwt")
def indexFastaGSE65525(infile, outfile):
    ''' build BWA index for merged fasta'''

    prefix = P.snip(outfile, ".1.ebwt")

    # build raw index
    statement = '''
    bowtie-build %(infile)s %(prefix)s >> %(outfile)s.log 2>&1
    '''
    P.run()


@follows(mkdir("GSE65525/transcriptome.dir"),
         indexFastaGSE65525)
@transform(processReadsGSE65525,
           regex("GSE65525/processed.dir/trimmed-(\S+)_UMI_(\S+).fastq.gz"),
           add_inputs(indexFastaGSE65525),
           r"GSE65525/transcriptome.dir/\1_UMI_\2.trans.bam")
def mapBowtieAgainstTranscriptomeGSE65525(infiles, outfile):
    ''' map reads using Bowtie against transcriptome

    bowtie parameterised according to Allon et al 2015 except
    random reporting of alignments where more than one "best" exist (-M 1):

    -n 1 number of mismatches allowed
    -l 15 seed length
    -e 300 maxmimum permitted sum of sequence qualities at all mismatched
           positions
    -M 1 if more than one "best" alignment exist, report one at random
    --best report in best to worst order
    --strata only report reads falling into the best stratum
    '''

    infile, reference = infiles
    job_threads = 2
    job_options = "-l mem_free=1.9G"
    bowtie_options = "-n1 -l 15 -e 300 -M 1 --best --strata"
    bowtie_index_dir = os.path.dirname(reference)
    genome = P.snip(os.path.basename(reference), ".1.ebwt")
    reffile = reference
    bowtie_threads = job_threads

    m = PipelineMapping.Bowtie(tool_options=bowtie_options,
                               remove_non_unique=0,
                               strip_sequence=0)

    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("GSE65525/dedup_unique.dir"),
         mkdir("GSE65525/dedup_percentile.dir"),
         mkdir("GSE65525/dedup_cluster.dir"),
         mkdir("GSE65525/dedup_adjacency.dir"),
         mkdir("GSE65525/dedup_directional_adjacency.dir"),
         mapBowtieAgainstTranscriptomeGSE65525)
@subdivide(mapBowtieAgainstTranscriptomeGSE65525,
           regex("GSE65525/transcriptome.dir/(\S+)_UMI_(\S+).trans.bam"),
           [r"GSE65525/dedup_unique.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE65525/dedup_percentile.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE65525/dedup_cluster.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE65525/dedup_adjacency.dir/\1_UMI_\2_deduped.trans.bam",
            r"GSE65525/dedup_directional_adjacency.dir/\1_UMI_\2_deduped.trans.bam"])
def dedupGSE65525(infile, outfiles):
    ''' perform deduping with various methods'''

    for outfile in outfiles:
        outfile_tmp = outfile.replace("_deduped.trans.bam", "_temp.trans.bam")
        outfile_snip = P.snip(outfile, ".bam")

        method = P.snip(os.path.basename(
            os.path.dirname(outfile).replace("dedup_", "")), ".dir")
        
        if method == "directional_adjacency":
            method = "directional-adjacency"

        if method == "adjacency":
            stats_cmd = ("--output-stats=%(outfile)s.stats "
                         "--further-stats" % locals())
        else:
            stats_cmd = ("--output-stats=%(outfile)s.stats " % locals())

        statement = '''umi_tools dedup
                   --per-contig --method=%(method)s %(stats_cmd)s
                   -I %(infile)s -S %(outfile_tmp)s -L %(outfile)s.log;
                   checkpoint;
                   samtools sort %(outfile_tmp)s %(outfile_snip)s;
                   checkpoint;
                   samtools index %(outfile)s;
                   rm -rf %(outfile_tmp)s'''
        P.run()


@collate([mapBowtieAgainstTranscriptomeGSE65525,
          dedupGSE65525],
         regex("GSE65525/(\S+).dir/(\S+)_UMI_(\S+).bam"),
         r"GSE65525/\1.dir/\2_UMI_\3.gene.counts.tsv")
def countGenesGSE65525(infile, outfile):
    ''' summarise counts per gene '''

    job_memory = "2G"
    PipelineScRNASeq.countAlignmentsPerGene(
        infile[0], outfile, mapq_threshold=10, submit=True)


@follows(mkdir("GSE65525/figures.dir"))
@collate(countGenesGSE65525,
         regex("GSE65525/(\S+).dir/(\S+)_UMI_(\S+).trans.gene.counts.tsv"),
         r"GSE65525/figures.dir/\1_\2_merged_gene_counts.tsv")
def mergeCountsGSE65525(infiles, outfile):
    ''' merge gene counts for each cell in each sample into a single
    table per sample'''

    job_memory = "4G"
    PipelineScRNASeq.getGeneCounts(infiles, outfile, submit=True)


@collate(dedupGSE65525,
         regex(r"GSE65525/(\S+)/(\S+)_UMI_\S+_deduped.trans.bam"),
         r"GSE65525/\1/\2_edit_distances.tsv")
def editDistanceDedupGSE65525(infiles, outfile):
    '''concatenate edit distance stats'''
    infiles = [infile + ".stats_edit_distance.tsv" for infile in infiles]

    job_memory = "1G"
    PipelineScRNASeq.summariseEditDistances(infiles, outfile, submit=True)


@follows(mkdir("GSE65525/figures.dir"))
@collate(editDistanceDedupGSE65525,
         regex("GSE65525/(\S+)/(\S+)_edit_distances.tsv"),
         [r"GSE65525/figures.dir/\2_merged_edit_distances.tsv",
          r"GSE65525/figures.dir/\2_edit_distances.png"])
def mergeAndPlotEditDistancesGSE65525(infiles, outfiles):
    '''merge edit distance stats for all deduplication methods and plot'''

    outfile, plot_out = outfiles

    job_memory = "1G"
    PipelineScRNASeq.mergeAndPlotEditDistances(
        infiles, outfile, plot_out, submit=True, job_memory=job_memory)


@collate(mergeCountsGSE65525,
         regex("GSE65525/figures.dir/(\S+)_SRR(\S+)_merged_gene_counts.tsv"),
         r"GSE65525/figures.dir/\1_merged_gene_counts.tsv")
def mergeGeneCountsPerDayGSE65525(infiles, outfile):
    ''' merge gene counts across timepoints, and add suffix to cell name '''

    job_memory = "12G"

    PipelineScRNASeq.mergeTimepointsGSE65525(
        infiles, outfile, submit=True, job_memory=job_memory)


@collate(mergeCountsGSE65525,
         regex("GSE65525/figures.dir/(\S+)_SRR(\S+)_merged_gene_counts.tsv"),
         r"GSE65525/figures.dir/SRR\2_cv.png")
def plotCVGSE65525(infiles, plotfile):
    ''' Calculate CV for each method and plot'''

    job_memory = "4G"

    normalise_method = "total-count"

    PipelineScRNASeq.plotCV(
        infiles, plotfile, normalise_method,
        submit=True, job_memory=job_memory)


@transform(mergeGeneCountsPerDayGSE65525,
           suffix(".tsv"),
           "_plots.log")
def clusterAndPlotGSE65525(infile, logfile):
    ''' '''
    pass


@merge(mergeGeneCountsPerDayGSE65525,
       "GSE65525/figures.dir/plots_heatmap.log")
def plotHeatmapsGSE65525(infiles, outfile):
    ''' plot heatmaps for each dedup method'''
    print infiles, outfile

    PipelineScRNASeq.plotHeatmapGSE65525(
        infiles, outfile, submit=True, job_memory="2G")


@transform(mergeGeneCountsPerDayGSE65525,
           regex("GSE65525/figures.dir/(\S+)_merged_gene_counts.tsv"),
           r"GSE65525/figures.dir/plots_PCA_\1_loadings.tsv")
def plotPCAGSE65525(infile, outfile):
    ''' plot PCA for each dedup method'''

    # big PCA -> ~ 5G memory required
    PipelineScRNASeq.plotPCAGSE65525(
        infile, outfile, submit=True, job_memory="5G")


@merge(plotPCAGSE65525,
       "GSE65525/figures.dir/PCA_variance.png")
def plotVariancePCAGSE65525(infiles, outfile):
    ''' plot varaince explained in first x PCs '''
    PipelineScRNASeq.plotVarianceGSE65525(
        infiles, outfile, submit=True)


@merge(plotPCAGSE65525,
       "GSE65525/figures.dir/PCA_loading.tsv")
def plotloadingsPCAGSE65525(infiles, outfile):
    ''' plot varaince explained in first x PCs '''
    PipelineScRNASeq.plotLoadingsGSE65525(
        infiles, outfile, submit=False)


@follows(mergeAndPlotEditDistancesGSE65525,
         mergeGeneCountsPerDayGSE65525,
         plotCVGSE65525,
         plotHeatmapsGSE65525,
         plotVariancePCAGSE65525,
         plotloadingsPCAGSE65525)
def GSE65525():
    pass


###############################################################################
# Section END
###############################################################################

# Combined edit distance plot
@mkdir("paper_figures.dir")
@merge((mergeAndPlotEditDistancesGSE65525,
        mergeAndPlotEditDistancesGSE53638),
       "paper_figures.dir/facetted_edit_distances.tsv")
def plotFacettedEditPlots(infiles, outfile):
    ''' combine the edit distance from all samples and plot the edit
    distances '''

    # we just want the first file (.tsv)
    infiles = [x[0] for x in infiles]

    PipelineScRNASeq.plotFacettedEditPlots(
        infiles, outfile, submit=False)


@follows(plotFacettedEditPlots)
def combinedPlots():
    pass


@follows(GSE65525,
         GSE53638,
         plotFacettedEditPlots)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
