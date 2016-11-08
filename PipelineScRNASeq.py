'''
PipelineScRNASeq.py - Utility functions for pipeline_scRNASeq.py
==============================================================

:Author: Toms Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Code
----

'''

import collections
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.Counts as Counts
import pysam
import pandas as pd
import numpy as np
import glob
import regex
import re
from CGATPipelines.Pipeline import cluster_runnable
import os
import CGAT.Fastq as Fastq
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
import pandas.rpy.common as com
import rpy2.robjects as robjects


@cluster_runnable
def extractUMIsAndFilterFastq(fastq_seq, fastq_UMI, barcodes):
    '''
    Paired end sequencing:
    read 1 - cell barcode (6bp) then UMI (10bp)
    read 2 - genomic sequence

    Need to extract UMI, append to read 2 read name and write out
    reads to single cell fastqs

    1. check sequence quality of barcode and UMI
    2. check barcode is an expected sequence (384 possible barcodes)
    3. Extract UMI from read 1 and append to read name for read pair
    4. write out to single cell fastq using barcode as identifier

    Filtering performed here mirrors filtering used by Soumillon et al 2014
    '''
    # TS - expected location of UMI and barcode is hard-coded in here...

    out_prefix = P.snip(fastq_UMI, "_1.fastq.gz")
    out_prefix_base = os.path.basename(out_prefix)

    fastq = IOTools.openFile(fastq_seq, "r")
    UMI_fastq = IOTools.openFile(fastq_UMI, "r")

    def FastqNext(infile):
        ''' get next fastq entry in file - without format checks'''
        line1 = infile.readline()
        if not line1:
            return None
        line2 = infile.readline()
        line3 = infile.readline()
        line4 = infile.readline()

        return Fastq.Record(line1[1:-1], line2[:-1], line4[:-1])

    record = FastqNext(fastq)
    UMI_record = FastqNext(UMI_fastq)

    def getReadName(record):
        return record.identifier.split(" ")[0]

    counts = collections.Counter()
    outfs = set()

    # use IOTools.FilePool class to open multiple file handles
    fastq_outfiles = IOTools.FilePool(
        output_pattern=out_prefix + "_UMI_%s.fastq.gz")

    while record is not None:
        # check both records are for the same fastq read
        assert getReadName(record) == getReadName(record), (
            "fastq read names do not match")

        UMI_record.format = "illumina-1.8"

        # check quality of barcode
        if not min(UMI_record.toPhred()[0:6]) >= 10:
            counts["barcode_quality_fail"] += 1

        else:
            # check quality of UMI
            if not min(UMI_record.toPhred()[6:16]) >= 30:
                counts["UMI_quality_fail"] += 1

            else:
                # check barcode matches supplied barcodes
                if not any(UMI_record.seq[0:6] == barcode
                           for barcode in barcodes):
                    counts["barcode_match_fail"] += 1

                else:
                    counts["keep"] += 1
                    cell = UMI_record.seq[0:6]
                    UMI = UMI_record.seq[6:16]
                    record.identifier = record.identifier.replace(" ", ":") + "_" + UMI
                    outfs.update((cell,))
                    fastq_outfiles.write(cell, "%s\n" % record)

        record = FastqNext(fastq)
        UMI_record = FastqNext(UMI_fastq)


@cluster_runnable
def countAlignmentsPerGene(infile, outfile, mapq_threshold=0):
    ''' count the number of reads aligning to each contig(gene) in bam '''
    # bash one-liner using cut|uniq -c would be
    # much quicker but the zero counts would be lost.

    counts = collections.Counter()
    insam = pysam.Samfile(infile, "rb")
    references = insam.references
    references_dict = {}

    for n, reference in enumerate(references):
        references_dict[n] = reference
        counts[reference] = 0

    inreads = insam.fetch()

    for read in inreads:
        if read.is_unmapped:
            continue
        if read.mapq < mapq_threshold:
            continue
        counts[references_dict[read.reference_id]] += 1

    with IOTools.openFile(outfile, "w") as f:
        for k, v in counts.most_common():
            f.write("%s\t%s\n" % (k, v))


@cluster_runnable
def summariseEditDistances(infiles, outfile):
    ''' tally counts for edit distances '''

    final_df = pd.DataFrame()
    for infile in infiles:
        temp_df = pd.read_table(infile, sep="\t", index_col="edit_distance")
        final_df = final_df.add(temp_df, fill_value=0)

    final_df = final_df.astype(int)
    final_df.to_csv(outfile, sep="\t")


@cluster_runnable
def mergeAndPlotEditDistances(infiles, outfile, plot_out):
    ''' merge the edit distance stats for each duplication method
    and plot cumulative sum'''

    plot_out2 = P.snip(plot_out, ".png") + "_unique_only.png"
    plot_out3 = P.snip(plot_out, ".png") + "_cluster_only.png"
    plot_out4 = P.snip(plot_out, ".png") + "_dir_only.png"

    final_df = pd.DataFrame()

    pre_null_df = pd.DataFrame()

    for infile in infiles:

        dedup_method = infile.split("/")[-2].replace(
            "dedup_", "").replace(".dir", "")

        tmp_df = pd.read_table(infile)
        tmp_df.columns = [x.replace("-", "_") for x in tmp_df.columns]
        tmp_post_df = tmp_df[[dedup_method, "%s_null" % dedup_method,
                              "edit_distance"]]
        tmp_post_df.columns = [x.replace(dedup_method, "post")
                               for x in tmp_post_df.columns]
        pre_null_df[dedup_method] = tmp_df["unique_null"]
        dedup_method = dedup_method.title()

        tmp_post_df['method'] = dedup_method

        final_df = final_df.append(tmp_post_df)

    final_df.to_csv(outfile, sep="\t", index=True)

    plot_edit_distances = R('''
    function(df){
    library(ggplot2)
    library(reshape)
    library(plyr)
    library(grid)

    df = df[df$edit_distance!="Single_UMI",]
    df$edit_distance = as.numeric(as.character(df$edit_distance))
    df_null = df[df$method=="Unique",]
    df_null$method="Null"
    df_null$post = df_null$post_null
    df = rbind(df, df_null)

    max_distance = max(df$edit_distance)
    col_range = topo.colors(max_distance, alpha = 0.9)

    df = df[df$edit_distance<max_distance,]
    df = df[order(df$edit_distance),]

    l_txt = element_text(size=20)

    df$method = factor(df$method,
                       levels=c("Unique", "Percentile", "Cluster",
                                "Adjacency", "Directional",
                                "Null"))

    a = aes(y=post, x=method)
    b = geom_bar(aes(fill=as.factor(edit_distance), colour="NA"),
                     stat="identity",
                     position="fill")
    x = xlab("")
    y = ylab("Fraction")
    s_f = scale_fill_manual(name="Edit distance", values=col_range)
    s_c = scale_colour_manual(guide="none", values="grey80")
    t= theme(axis.text.y=l_txt,
             axis.title.x=l_txt, axis.title.y=l_txt,
             legend.text=l_txt, strip.text=l_txt,
             legend.title=l_txt,
             panel.margin = unit(1, "lines"),
             legend.position="right",
             legend.key.size=unit(2, "line"),
             aspect.ratio=1,
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

    x_no_angle = theme(axis.text.x=element_text(size=20))
    x_angle = theme(axis.text.x=element_text(size=20, angle=90,
                                             hjust=1, vjust=0.5))

    p = ggplot(df) + a + b + x + y + s_f + theme_bw() + t + x_angle + s_c

    ggsave("%(plot_out)s", width=10, height=10)

    p2 = ggplot(df[(df$method=="Unique") | (df$method=="Null"),]) +
         a + b + x + y + theme_bw() + s_f + t + x_no_angle + s_c

    ggsave("%(plot_out2)s", width=10, height=8)

    p3 = ggplot(df[(df$method=="Cluster") | (df$method=="Null"),]) +
         a + b + x + y + theme_bw() + s_f + t + x_no_angle + s_c

    ggsave("%(plot_out3)s", width=10, height=8)

    p4 = ggplot(df[(df$method=="Directional") | (df$method=="Null"),]) +
         a + b + y + x + theme_bw() + s_f + t + x_no_angle + s_c

    ggsave("%(plot_out4)s", width=10, height=8)
    }
    ''' % locals())

    final_df.index = range(0, len(final_df.index), 1)
    r_df = com.convert_to_r_dataframe(final_df)
    plot_edit_distances(r_df)


@cluster_runnable
def getGeneCounts(infiles, outfile):

    df = pd.DataFrame()

    for infile in infiles:
        temp_df = pd.io.parsers.read_csv(
            infile, sep="\t", index_col=0, names=["counts"])

        cell = re.sub(".*_UMI_", "", infile).replace(
            "_deduped", "").replace(
            ".trans.gene.counts.tsv", "")
        temp_df.index = temp_df.index
        temp_df.sort_index(inplace=True)
        temp_df.index.names = ['gene']
        df[cell] = temp_df["counts"]

    df.to_csv(outfile, sep="\t", index=True)


@cluster_runnable
def makeHeatmapsAndPCA(infiles,
                       outfile,
                       min_counts_per_sample=1000,
                       max_counts_per_sample=100000,
                       top_genes_heatmap=100, top_genes_pca=2000):
    '''
    1. Concatenates the counts from day 0 and day 14 and subset to samples
       with counts within set range.
    2. Heatmaps plotted for each dedup method using these samples
       and top 100 genes
    3. PCA plotted for each dedup method using these samples and top 2000 genes
    '''

    def getDay(infile):
        ''' return the day from the infile sample id'''
        if "SRR1058003" in infile:
            return "_day0"

        elif "SRR1058023" in infile:
            return "_day14"

        else:
            return None

    methods = set([re.sub("_SRR.*", "", os.path.basename(x).replace(
        "dedup_", "")) for x in infiles])

    uniq_final_df = pd.DataFrame()

    with IOTools.openFile(outfile, "w") as outf:

        # start with just the uniq
        for infile in [x for x in infiles if "unique" in x]:

            day = getDay(infile)

            # ignore sample id SRR1058003 / 1058023
            if not day:
                continue

            tmp_df = pd.io.parsers.read_csv(infile, sep="\t", index_col=0)
            tmp_df.columns = [x + day for x in tmp_df.columns]

            uniq_final_df = pd.concat([uniq_final_df, tmp_df], axis=1)

        uniq_genes_only = uniq_final_df.drop(
            [x for x in uniq_final_df.index if "ERCC" in str(x)])
        uniq_genes_only = uniq_genes_only.drop("chrM")

        # filter out low abundance and high abundance samples
        size_factors = uniq_genes_only.sum(axis=0)

        keep = [(float(x) > min_counts_per_sample and
                 float(x) < max_counts_per_sample) for x in size_factors]

        uniq_filtered = uniq_genes_only.iloc[:, keep]
        cells = uniq_filtered.columns
        outf.write("number of cells: %i\n" % len(cells.tolist()))

        # need to recalculate as some samples have been removed
        size_factors = uniq_filtered.sum(axis=0)
        r_size_factors = robjects.Vector(size_factors.tolist())

        # normalise
        uniq_normed = uniq_filtered * 1000000.0 / size_factors

        # subset to option(**top_genes**) most highly expressed
        sum_counts_per_row = uniq_normed.sum(1)

        countsLog = Counts.Counts(uniq_normed)
        countsLog.normalise(method="total-column")
        countsLog.log()
        log_df = countsLog.table

        countsZ = Counts.Counts(uniq_normed)
        countsZ.zNormalise()
        z_df = countsZ.table

        z_df['order'] = sum_counts_per_row
        log_df['order'] = sum_counts_per_row

        z_df = z_df.sort(columns="order", ascending=False)
        log_df = log_df.sort(columns="order", ascending=False)

        z_df_pca = z_df[0:top_genes_pca]
        z_df_heatmap = z_df[0:top_genes_heatmap]
        log_df_heatmap = log_df[0:top_genes_heatmap]

        z_df_pca_filtered = z_df_pca.drop("order", 1)
        z_df_heatmap_filtered = z_df_heatmap.drop("order", 1)
        log_df_heatmap_filtered = log_df_heatmap.drop("order", 1)

        genes_pca = z_df_pca_filtered.index
        genes_heatmap = z_df_heatmap_filtered.index
        genes_heatmap_log = z_df_heatmap_filtered.index
        outf.write("number of pca genes: %i\n" % len(genes_pca.tolist()))
        outf.write("number of heatmap genes: %i\n" % len(genes_heatmap.tolist()))

        plot_heatmap = R('''
        function(df, plot_outfile, table_outfile){
        options(expressions = 10000)
        library(Biobase)
        library(RColorBrewer)
        library(gplots)

        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

        col_cols = sapply(colnames(df),
          FUN=function(x) ifelse(grepl(".*day0", x), "chartreuse4", "chocolate2"))

        d = as.dendrogram(hclust(as.dist(1-cor(df,
                                       method = "spearman")),
                         method = "ward.D2"))

        # write out clusters
        df_clusters = data.frame(
            "cluster" = dendextend:::cutree.dendrogram(d, k=2))
        df_clusters$sample_group = gsub("[A-Z]+_","", rownames(df_clusters))
        df_clusters$sample_group <- gsub("day0", "mES Cells",
                                         df_clusters$sample_group)
        cluster_table = table(df_clusters$cluster, by=df_clusters$sample_group)

        write.table(cluster_table, table_outfile,
                    quote=FALSE, sep="\t", row.names=FALSE)

        png(plot_outfile, width=25, height=25, unit="cm", res=400)

        heatmap.2(as.matrix(df),
                  col = hmcol, scale="none", trace="none",
                  margin=c(18, 10), keysize=1, cexCol=2,
                  dendrogram="column",
                  Colv = d,
                  ColSideColors=col_cols,
                  hclustfun = function(x) hclust(x, method = 'average'),
                  distfun = function(x) as.dist(1 - cor(t(x),
                                                method="spearman")),
                  labRow=F, labCol=F, key=F)

        legend(0, 1,
        legend = unique(sapply(colnames(df),
                 FUN=function(x) ifelse(grepl(".*_day0", x),
                                        "Day 0", "Day 14"))),
        col = unique(col_cols),
        lty=1, lwd=0, pch=15, cex=1, pt.cex = 3, bty="n")

        dev.off()

        }''')

        plot_PCA = R('''
        function(df,
                 r_size_factors,
                 plot_base){

        library(ggplot2)
        library(grid)

        m_text = element_text(size=15)
        s_text = element_text(size=10)

        gene_pca <- prcomp(t(df), center = TRUE)
        variance = gene_pca$sdev^2
        variance_explained = round(variance/sum(variance), 5)

        variance_df = data.frame("Variance_explained" = variance_explained,
                                 "PC" = seq(1, length(variance)))

        write.table(variance_df, paste0(plot_base, "_variance.tsv"),
                    sep="\t", quote=FALSE, row.names=FALSE)

        p_variance = ggplot(variance_df, aes(x=PC, y=Variance_explained))+
        geom_point()+
        geom_line()+
        theme_classic()+
        ylab("Variance explained (%%)")+
        theme_bw() +
        theme(axis.text.x = m_text,
            axis.title.y = m_text,
            axis.title.x = m_text,
            axis.text.y = m_text,
            aspect.ratio=1)

        ggsave(paste0(plot_base, "_pca_variance.png"), width=5, height=5)

        PCs_df = data.frame(gene_pca$x)
        PCs_df$id_1 = sapply(strsplit(rownames(PCs_df), "_"), "[", 1)
        PCs_df$id_2 = sapply(strsplit(rownames(PCs_df), "_"), "[", 2)
        PCs_df$id_expression = log(as.numeric(r_size_factors),10)

        write.table(PCs_df, paste0(plot_base, "_eigenvectors.tsv"),
                    sep="\t", row.names=FALSE)

        p = geom_point(size=2, aes(shape=id_2, colour=id_2))
        s_c = scale_colour_discrete(name="")
        s_s = scale_shape_discrete(name="")
        t = theme(axis.text.x = s_text, axis.text.y = s_text,
                  title = m_text, legend.text = m_text,
                  legend.title = m_text, legend.key.size = unit(7, "mm"),
                  aspect.ratio=1)

        p_pca = ggplot(PCs_df, aes(x=PC1, y=PC2)) +
        xlab(paste0('PC1 (Variance explained = ' ,
                    round(100 * variance_explained[1], 1), '%)')) +
        ylab(paste0('PC2 (Variance explained = ' ,
                     round(100 * variance_explained[2], 1), '%)')) +
        p + theme_bw() + s_c + s_s + t

        ggsave(paste0(plot_base, "_pc1_pc2.png"), width=5, height=5)

        p_pca = ggplot(PCs_df, aes(x=PC3, y=PC4)) +
        xlab(paste0('PC4 (Variance explained = ' ,
                    round(100 * variance_explained[3], 1), '%)')) +
        ylab(paste0('PC4 (Variance explained = ' ,
                     round(100 * variance_explained[4], 1), '%)')) +
        p + theme_bw() + s_c + s_s + t

        ggsave(paste0(plot_base, "_pc3_pc4.png"), width=5, height=5)

        p_pca = ggplot(PCs_df, aes(x=PC1, y=PC2)) +
        geom_point(size=2, aes(shape=id_2, colour=id_expression)) +
        scale_shape_manual(name=guide_legend(title='Sample'), values=c(8,20)) +
        scale_colour_continuous(name=guide_legend(title='Counts')) +
        xlab(paste0('PC1 (Variance explained = ' ,
                    round(100 * variance_explained[1], 1), '%)')) +
        ylab(paste0('PC2 (Variance explained = ' ,
                     round(100 * variance_explained[2], 1), '%)')) +
        theme_bw() + t

        ggsave(paste0(plot_base, "_pc1_pc2_expression.png"), width=5, height=5)

        p_pca = ggplot(PCs_df, aes(x=PC3, y=PC4)) +
        geom_point(size=2, aes(shape=id_2, colour=id_expression)) +
        scale_shape_manual(name=guide_legend(title='Sample'), values=c(8,20)) +
        scale_colour_continuous(name=guide_legend(title='Counts')) +
        xlab(paste0('PC3 (Variance explained = ' ,
                    round(100 * variance_explained[1], 1), '%)')) +
        ylab(paste0('PC4 (Variance explained = ' ,
                     round(100 * variance_explained[2], 1), '%)')) +
        theme_bw() + t

        ggsave(paste0(plot_base, "_pc3_pc4_expression.png"), width=7, height=5)

        }''')

        plot_outfile = P.snip(outfile, ".log") + "_uniq.png"
        cluster_table = P.snip(outfile, ".log") + "_clusters.tsv"

        plot_base = os.path.join(os.path.dirname(outfile), "uniq")
        outf.write("%s\n" % plot_base)

        r_uniq_pca_df = com.convert_to_r_dataframe(z_df_pca_filtered)
        r_uniq_heatmap_df = com.convert_to_r_dataframe(log_df_heatmap_filtered)

        log_df_heatmap_filtered.to_csv(outfile+"4", sep="\t", index=True)
        plot_heatmap(r_uniq_heatmap_df, plot_outfile, cluster_table)
        plot_PCA(r_uniq_pca_df,
                 r_size_factors,
                 plot_base)

        outf.write("plotting heatmap for uniq deduping\n")
        outf.write("plotting PCAs for uniq deduping\n")

        # now repeat the heatmap and pca plotting for the other methods
        outf.write("methods:\n%s\n" % "\n".join(map(str, methods)))
        methods.remove("unique")
        methods.remove("transcriptome")

        for method in methods:

            final_df = pd.DataFrame()

            for infile in infiles:

                if method == "transcriptome":
                    pattern = method
                else:
                    pattern = "dedup_%s" % method

                if pattern in infile:
                    outf.write("%s\n" % infile)
                    day = getDay(infile)

                    # ignore other sample ids
                    if not day:
                        continue

                    tmp_df = pd.io.parsers.read_csv(infile, sep="\t", index_col=0)
                    tmp_df.columns = [x + day for x in tmp_df.columns]

                    final_df = pd.concat([final_df, tmp_df], axis=1)

            outf.write("%s\n" % "\t".join(map(str, final_df.shape)))
            outf.write("number of cells: %i\n" % len(final_df.index))
            # subset to chosen cells from uniquq dedup and
            # extract the size factors for the selected cells
            filtered = final_df.loc[:, cells.tolist()]
            outf.write("%s\n" % "\t".join(map(str, filtered.shape)))

            size_factors = filtered.sum(axis=0)
            r_size_factors = robjects.Vector(size_factors.tolist())

            normed = filtered * 1000000.0 / size_factors

            countsLog = Counts.Counts(normed)
            countsLog.normalise(method="total-column")
            countsLog.log(inplace=True)

            countsZ = Counts.Counts(normed)
            countsZ = countsZ.zNormalise(inplace=False)

            z_df = countsZ.table
            log_df = countsLog.table

            # subset to genes from unique dedup
            normed_pca_filtered = z_df.loc[genes_pca, :]
            normed_heatmap_filtered = log_df.loc[genes_heatmap, :]
            outf.write("%s\n" % "\t".join(map(str, normed_pca_filtered.shape)))
            outf.write("%s\n" % "\t".join(map(str, normed_heatmap_filtered.shape)))

            plot_outfile = P.snip(outfile, ".log") + "_%s.png" % method
            cluster_table = P.snip(outfile, ".log") + "_clusters_%s.tsv" % method
            plot_base = os.path.join(os.path.dirname(outfile), method)

            normed_pca_filtered.to_csv(
                outfile.replace(".log", "%s_pca_df.tsv" % method), sep="\t")
            normed_heatmap_filtered.to_csv(
                outfile.replace(".log", "%s_heatmap_df.tsv" % method), sep="\t")

            r_df_pca = pandas2ri.py2ri(normed_pca_filtered)
            r_df_heatmap = pandas2ri.py2ri(normed_heatmap_filtered)

            outf.write("plotting heatmap for %s deduping\n" % method)
            plot_heatmap(r_df_heatmap, plot_outfile, cluster_table)

            outf.write("plotting PCAs for %s deduping\n" % method)
            plot_PCA(r_df_pca,
                     r_size_factors,
                     plot_base)


@cluster_runnable
def extractUMIsAndFilterFastqGSE65525(fastq_UMI, fastq_seq,
                                      barcodes1_infile, barcodes2_infile,
                                      cell_barcode_count):
    '''Paired end sequencing:
    read 1 - 51bp: cell barcode1 (8-12bp) then adapter sequence (22bp),
             cell barcode2 (8bp) then UMI (6bp), then Ts
    read 2 - genomic sequence

    Need to extract UMI, append to read 2 read name and write out
    reads to single cell fastqs for each cell

    1. Check barcode is an expected sequence (384^2 possible barcodes)
    2. Extract UMI and cell barcode from read 1 and append to read
       name for read pair2
    3. Write out to temporary single end fastq and generate
       frequency table of barcodes
    4. Identify n most abundant barcodes where n is the the of cells
       and re-parse single-end fastq to retain only these barcodes
    5. Strip cell barcode from fastq identifier and write out to
       single cell fastq using concatenated barcodes as fastq name

    Filtering performed here mirrors filtering used by Klein et al
    2015

    Ideally, we would parse the pair of fastqs twice. First to
    generate the frequency table and then again to identify all cell
    barcodes within a hamming distance of 2 from the n most abundnant
    barcodes. However, the identification of cell barcode near
    mismatches is extremely time consuming and some barcodes cannot be
    unambiguously resolved at a Hamming distance of 2. Hence, only
    perfect cell barcode matches are retained.

    Since we have already identified all the candidate perfect matches
    on the first parse, we can save these out to a temporary single
    end fastq and parse this to identify the perfect matches to the n
    most abundant cell barcodes
    '''

    def reverseComp(seq):
        comp = {"A": "T",
                "T": "A",
                "C": "G",
                "G": "C",
                "N": "N"}

        return "".join([comp[x] for x in seq[::-1]])

    barcode_set1 = set()
    barcode_set2 = set()

    with IOTools.openFile(barcodes1_infile) as inf:
        for line in inf:
            barcode_set1.update((reverseComp(line.strip()),))

    with IOTools.openFile(barcodes2_infile) as inf:
        for line in inf:
            barcode_set2.update((reverseComp(line.strip()),))

    out_prefix = P.snip(fastq_UMI, "_1.fastq.gz")
    #out_prefix_base = os.path.basename(out_prefix)

    fastq = IOTools.openFile(fastq_seq, "r")
    UMI_fastq = IOTools.openFile(fastq_UMI, "r")

    def FastqNext(infile):
        ''' get next fastq entry in file - without format checks'''
        line1 = infile.readline()
        if not line1:
            return None
        line2 = infile.readline()
        line3 = infile.readline()
        line4 = infile.readline()

        return Fastq.Record(line1[1:-1], line2[:-1], line4[:-1])

    def getReadName(record):
        return record.identifier.split(" ")[0]

    counts = collections.Counter()
    outfs = collections.Counter()

    outfile_tmp = out_prefix + "_tmp.fastq"
    fastq_tmp_outfile = IOTools.openFile(outfile_tmp, "w")

    n = 0
    with IOTools.openFile(out_prefix + "_log.tsv", "w") as outf:
        # should really add an assert here to check for empty files
        record = 1
        while record is not None:

            record = FastqNext(fastq)
            UMI_record = FastqNext(UMI_fastq)

            # check for end of file
            if record is None or UMI_record is None:
                if UMI_record is None:
                    outf.write("reached end of UMI fastq\n")
                if record is None:
                    outf.write("reached end of paired end fastq\n")
                break

            if n % 100000 == 0:
                outf.write("read through %i fastq records. %i records retained"
                           "\n" % (n, counts["keep"]))
                for reason, count in counts.most_common():
                    outf.write("%s\t%i\n" % (reason, count))

            n += 1

            # check both records are for the same fastq read
            assert getReadName(UMI_record) == getReadName(record), (
                "fastq read names do not match")

            UMI_record.format = "illumina-1.8"

            # match must start from at least 9 bp into fastq sequence
            # allow two mismatches
            pos = regex.search("(GAGTGATTGCTTGTGACGCCTT){s<=2}",
                               UMI_record.seq[8:])
            
            #pos = re.search("GAGTGATTGCTTGTGACGCCTT", UMI_record.seq)

            # Skip if:
            # - no match
            # - matched string is not the right length
            # - match starts too far from 5' end
            if not pos:
                counts["no_adapter_sequence"] += 1
                continue

            if len(pos.captures()[0]) != 22:
                counts["adapter_sequence_wrong_length"] += 1
                continue
                
            start = pos.start() + 8
            if start > 12:
                counts["no_adapter_sequence_match_in_correct_location"] += 1
                continue

            barcode1 = UMI_record.seq[0:start]
            barcode2 = UMI_record.seq[start+22:start+30]
            UMI = UMI_record.seq[start+30:start+36]

            # check barcode and UMI lengths. This is mainly to exlude
            # read lengths which are too short to encode the UMI but
            # also a back up check for barcode 1
            if len(barcode1) < 8:
                counts["barcode1_too_short"] += 1
                continue

            if len(barcode1) > 12:
                counts["barcode1_too_long"] += 1
                continue

            if len(barcode2) != 8:
                counts["barcode2_wrong_length"] += 1
                continue

            if len(UMI) < 6:
                counts["UMI_too_short"] += 1
                continue

            # check barcodes are found in lists
            if str(barcode1) not in barcode_set1:
                counts["barcode1_mismatch"] += 1
                continue

            if str(barcode2) not in barcode_set2:
                counts["barcode2_mismatch"] += 1
                continue

            counts["kept"] += 1
            cell = barcode1 + barcode2

            # for now, just write out all the fastqs which pass filters
            # we'll reparse through this fastq later and "correct" the cell barcode
            record.identifier = "_".join((
                record.identifier.replace(" ", ":"), UMI,  cell))

            outfs[cell] += 1
            fastq_tmp_outfile.write("%s\n" % record)

        for reason, count in counts.most_common():
            outf.write("%s\t%i\n" % (reason, count))

    fastq_tmp_outfile.close()
    fastq.close()
    UMI_fastq.close()

    with IOTools.openFile(out_prefix + "_barcode_counts.tsv", "w") as outf:
        for barcode, count in outfs.most_common():
            outf.write("%s\t%i\n" % (barcode, count))

    with IOTools.openFile(out_prefix + "_barcode1.tsv", "w") as outf:
        outf.write("barcode set 1\n")
        for x in barcode_set1:
            outf.write("%s\n" % "\t".join((x, "length: ", str(len(x)))))

    with IOTools.openFile(out_prefix + "_barcode2.tsv", "w") as outf:
        outf.write("barcode set 2\n")
        for x in barcode_set2:
            outf.write("%s\n" % "\t".join((x, "length: ", str(len(x)))))

    # find the first n cell barcodes these represent the selected cells
    cell_barcode_set = set()
    n = 1
    for barcode, count in outfs.most_common():
        if n <= cell_barcode_count:
            n += 1
            cell_barcode_set.update((barcode,))
        else:
            break

    fastq = IOTools.openFile(outfile_tmp, "r")

    # use IOTools.FilePool class to open multiple file handles
    fastq_outfiles = IOTools.FilePool(
        output_pattern=out_prefix + "_UMI_%s.fastq.gz")

    # re-initate objects
    counts = collections.Counter()
    outfs = collections.Counter()

    counts["total"] = 0
    with IOTools.openFile(out_prefix + "_log2.tsv", "w") as outf:
        # should really add an assert here to check for empty files
        record = 1
        while record is not None:

            record = FastqNext(fastq)

            # check for end of file
            if record is None:
                outf.write("reached end of fastq\n")
                break

            if counts["total"] % 100000 == 0:
                outf.write("read through %i fastq records. %i records retained"
                           "\n" % (counts["total"], counts["kept"]))
                for reason, count in counts.most_common():
                    outf.write("%s\t%i\n" % (reason, count))

            counts["total"] += 1

            record.format = "illumina-1.8"

            cell = record.identifier.split("_")[-1]

            # check if cell barcode is in set
            if cell not in cell_barcode_set:
                counts["cell_barcode_mismatch"] += 1
                continue

            else:
                counts["kept"] += 1

            # append the UMI onto read2 and write out
            #record.identifier = "_".join((
            #    record.identifier.replace(" ", ":"), UMI))
            record.identifier = "_".join((
                record.identifier.split("_")[0:2]))

            outfs[cell] += 1
            fastq_outfiles.write(cell, "%s\n" % record)

        for reason, count in counts.most_common():
            outf.write("%s\t%i\n" % (reason, count))

    with IOTools.openFile(out_prefix + "_barcode_counts2.tsv", "w") as outf:
        for barcode, count in outfs.most_common():
            outf.write("%s\t%i\n" % (barcode, count))

    os.unlink(outfile_tmp)


@cluster_runnable
def mergeTimepointsGSE65525(infiles, outfile):
    ''' merge gene counts from different timepoints and write out '''

    df = pd.DataFrame()

    for infile in infiles:
        if "SRR1784310" in infile:
            day = 0
        elif "SRR1784313" in infile:
            day = 2
        elif "SRR1784314" in infile:
            day = 4
        elif "SRR1784315" in infile:
            day = 7

        tmp_df = pd.read_csv(IOTools.openFile(infile, "r"),
                             sep="\t", index_col=0)

        tmp_df.columns = ["%s_day%i" % (x, day) for x in tmp_df.columns]

        df = pd.concat([df, tmp_df], axis=1)

    df.to_csv(outfile, sep="\t", index=True)


@cluster_runnable
def plotCV(infiles, plotfile, normalise_method):
    ''' calculate CVs and plot, split by method'''

    plotfile2 = plotfile.replace(".png", "_difference.png")
    plotfile3 = plotfile.replace(".png", "_no_unique.png")
    plotfile4 = plotfile.replace(".png", "_vs_mean.png")
    plotfile5 = plotfile.replace(".png", "_difference_hist.png")

    cv_df = pd.DataFrame()

    for infile in infiles:
        method = re.sub("_SRR.*", "", os.path.basename(infile)).replace(
            "dedup_", "")
        method = method.title()

        if method == "Transcriptome":
            method = "None"

        tmp_df = Counts.Counts(pd.read_csv(IOTools.openFile(infile, "r"),
                                           sep="\t", index_col=0))

        tmp_df.removeObservationsFreq(1)
        tmp_df.normalise(method=normalise_method)

        tmp_cv_df = pd.DataFrame({
            "Method": method,
            "mean": tmp_df.table.apply(np.mean, axis=1),
            "cv": tmp_df.table.apply(
                func=lambda(x): np.std(x)/np.mean(x), axis=1)})

        cv_df = pd.concat([cv_df, tmp_cv_df], axis=0)

    # need to check all genes have CVs for all methods
    # it's possible for some genes to have no CV value for network-based methods
    # if all reads fail the mapq threshold --> all zero counts --> no CV!
    # we can do this by operating on the index = gene names

    n_methods = len(set(cv_df['Method']))
    gene_counts = collections.Counter(cv_df.index)
    retain_genes = [x for x in gene_counts.keys() if gene_counts[x] == n_methods]
    cv_df = cv_df.ix[retain_genes]

    table_outfile = plotfile.replace(".png", ".tsv")
    cv_df.to_csv(table_outfile, sep="\t", index=True)

    plotCV = R('''
    function(){

    library(ggplot2)

    df = read.table("%(table_outfile)s", sep="\t", header=T)

    m_txt = element_text(size=18)

    min_x = log(min(df$mean), 10)
    max_x = log(max(df$mean), 10)
    min_y = log(min(df$cv), 10)
    max_y = log(max(df$cv), 10)

    df$Method = factor(df$Method, levels=c("None", "Unique", "Percentile",
                                           "Cluster", "Adjacency",
                                           "Directional"))

    t = theme(axis.text.x = m_txt,
              axis.text.y = m_txt,
              axis.title.x = m_txt,
              axis.title.y = m_txt,
              legend.position = "none",
              strip.text = m_txt,
              aspect.ratio=1)


    p = ggplot(df, aes(x = log(mean, 10) , y = log(cv, 10))) +
    geom_point(size=1, alpha=0.2) +
    geom_line(colour="red", aes(x=log(mean, 10),
                                y = log(sqrt(mean)/mean, 10))) +
    facet_wrap(~Method,nrow=2) +
    theme_bw() + t +
    xlim(min_x, max_x) +
    ylim(min_y, max_y) +
    xlab("Mean expression (log10)") +
    ylab("CV (log10)")

    ggsave("%(plotfile)s", width=15, height=10)

    df = df[order(df$Method),]
    df["unique_cv"] = df[df['Method']=="Unique", "cv"]

    p = ggplot(df, aes(log(unique_cv,10), log(cv,10))) +
    geom_point(size=2, aes(colour=Method)) +
    geom_line(aes(x=log(unique_cv, 10), y=log(unique_cv,10)),
              colour="black", size=0.25) +
    facet_wrap(~Method,nrow=2) + theme_bw() + t +
    xlab("'Unique' CV (log10)") +
    ylab("CV (log10)")

    ggsave("%(plotfile2)s", width=15, height=10)

    methods_keep = c("Percentile", "Cluster",
                     "Adjacency", "Directional")

    df = df[df$Method %%in%% methods_keep,]

    p = ggplot(df, aes(log(unique_cv,10), log(cv,10), colour=Method)) +
    geom_point(alpha=0.6) +
    facet_wrap(~Method,nrow=2) +
    theme_bw() + t +
    geom_line(size=0.5, colour="black", aes(log(unique_cv,10),log(unique_cv,10))) +
    xlab("'Unique' CV) (log10)") +
    ylab("CV (log10)")

    ggsave("%(plotfile3)s", width=10, height=10)

    min_y_diff = min(log(df$cv-df$unique))
    max_y_diff = max(log(df$cv-df$unique))

    df$binned_mean = .bincode(log(df$mean,10), c(seq(-4,2,1), 10))
    p = ggplot(df, aes(x = binned_mean, y = cv-unique_cv)) +
    geom_boxplot() + #(size=1, alpha=0.2) +
    facet_wrap(~Method,nrow=2) +
    theme_bw() + t +
    xlab("Mean expression (log10)") +
    ylab("CV difference (relative to 'Unique' CV)")

    ggsave("%(plotfile4)s", width=10, height=10)

    t2 = theme(axis.text.x = m_txt,
              axis.text.y = m_txt,
              axis.title.x = m_txt,
              axis.title.y = m_txt,
              strip.text = m_txt,
              aspect.ratio=1,
              legend.position="none")

    p = ggplot(df, aes((cv - unique_cv))) +
    geom_histogram(aes(fill=Method), breaks=seq(-5.1,5.1,0.2)) +
    xlim(-2.5, 2.5) +
    facet_wrap(~Method, nrow=2) +
    theme_bw() + t2 +
    xlab("CV difference (relative to Unique CV)") +
    ylab("Count")

    ggsave("%(plotfile5)s", width=7, height=7)

    }''' % locals())

    plotCV()


@cluster_runnable
def clusterAndPlot(infile, outfile_prefix):
    ''' cluster counts and plot uses hierachical clustering
    (spearman's rank, average linkage) on top 1000 most highly
    expressed genes. Plot heatmap with dendograms'''

    counts = Counts.Counts(pd.read_csv(IOTools.openFile(infile, "r"),
                                       sep="\t", index_col=0))

    # note: can't use DESeq size factors as the median of many columns is zero
    # even after removing rows with just zeros! use naive normalisation instead
    counts.removeObservationsFreq(min_counts_per_row=1)
    counts.normalise("million-counts")
    counts.log(inplace=True)


@cluster_runnable
def plotHeatmapGSE65525(infiles,
                        outfile,
                        percentile=90):
    ''' '''

    plot_heatmap = R('''
        function(df, plot_outfile, table_outfile){

        options(expressions = 10000)

        library(Biobase)
        library(RColorBrewer)
        library(gplots)
        library("dendextend")

        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

        col_cols = sapply(colnames(df),
        FUN=function(x) ifelse(grepl(".*day0", x), "grey50",
                            ifelse(grepl(".*day2", x), "#E69F00"
                                ,ifelse(grepl(".*day4", x), "#56B4E9",
                                        "#F0E442"))))

        d = as.dendrogram(hclust(as.dist(1-cor(log(df+0.1),
                                       method = "spearman")),
                         method = "ward.D2"))

        colorCodes <- c(day0="grey50", day2="#E69F00",
                        day4="#56B4E9", day7="#F0E442")

        groupCodes <- gsub("[A-Z]+_","",colnames(df))
        labels_colors(d) <- colorCodes[groupCodes][order.dendrogram(d)]

        # write out clusters
        df_clusters = data.frame(
            "cluster" = dendextend:::cutree.dendrogram(d, k=4))
        df_clusters$sample_group = gsub("[A-Z]+_","", rownames(df_clusters))
        df_clusters$sample_group <- gsub("day0", "mES Cells",
                                         df_clusters$sample_group)
        cluster_table = table(df_clusters$cluster, by=df_clusters$sample_group)

        write.table(cluster_table, table_outfile,
                    quote=FALSE, sep="\t", row.names=FALSE)

        png(gsub(".png", "_dendo.png", plot_outfile), width=25, height=25,
            unit="cm", res=400)
        plot(d)
        dev.off()

        png(plot_outfile, width=25, height=25, unit="cm", res=400)

        heatmap.2(as.matrix(log(df+0.1)),
          col = hmcol, scale="none", trace="none",
          margin=c(18, 10), keysize=1, cexCol=2,
          dendrogram="column",
          Colv = d,
          ColSideColors=col_cols,
          hclustfun = function(x) hclust(x, method = 'average'),
          distfun = function(x) as.dist(1 - cor(t(x),
                                        method="spearman")),
          labRow=F, labCol=F, key=F)

        legend(0, 1,
        legend = unique(sapply(colnames(df),
            FUN=function(x) ifelse(grepl(".*day0", x), "mES Cells",
                                ifelse(grepl(".*day2", x), "Day 2",
                                       ifelse(grepl(".*day4", x), "Day 4", "Day 7"))))),
        col = unique(col_cols),
        lty=1, lwd=0, pch=15, cex=1, pt.cex=3, bty="n")

        dev.off()

        }''')

    with IOTools.openFile(outfile, "w") as outf:

        # start with just the uniq
        for infile in [x for x in infiles if "unique" in x]:

            uniq_df = pd.read_table(infile, sep="\t", index_col=0)

            counts = Counts.Counts(uniq_df)
            counts.normalise("total-count")

            df = counts.table
            counts.removeObservationsPerc(percentile_rowsums=percentile)
            genes = counts.table.index

            top_df = com.convert_to_r_dataframe(counts.table)

            plot_outfile = P.snip(outfile, ".log") + "_heatmap_uniq.png"
            plot_table = P.snip(outfile, ".log") + "_dendogram_clusters_uniq.tsv"
            uniq_outfile = P.snip(outfile, ".log") + "_table_uniq.tsv"

            plot_heatmap(top_df, plot_outfile, plot_table)
            counts.table.to_csv(uniq_outfile, sep="\t")
            outf.write("made unique heatmap\n %s" % plot_outfile)

        outf.write("%s\n" % "\n".join(infiles))
        # now do the rest of the methods using the genes in unique
        for infile in [x for x in infiles if "unique" not in x]:
            method = os.path.basename(infile).replace("dedup_", "").replace(
                "_merged_gene_counts.tsv", "")

            df = pd.io.parsers.read_csv(infile, sep="\t", index_col=0)

            counts = Counts.Counts(df)
            counts.normalise("total-count")
            df = counts.table
            counts.removeObservationsPerc(percentile_rowsums=percentile)
            top_df = com.convert_to_r_dataframe(counts.table)

            plot_outfile = P.snip(outfile, ".log") + "_heatmap_%s.png" % method
            plot_table = P.snip(outfile, ".log") + "_dendogram_clusters_%s.tsv" % method
            method_outfile = P.snip(outfile, ".log") + "_table_%s.tsv" % method

            counts.table.to_csv(method_outfile, sep="\t")
            plot_heatmap(top_df, plot_outfile, plot_table)
            outf.write("made %s heatmap\n" % method)


@cluster_runnable
def plotPCAGSE65525(infile,
                    pca_loadings):
    ''' '''

    plot_PCA = R('''

    function(df, variance_table, variance_outfile, vector_outfile,
             pca_outfile1, pca_outfile2, pca_loadings){

    library(ggplot2)

    gene_pca <- prcomp(t(df), center = TRUE)

    m_text = element_text(size=15)
    s_text = element_text(size=10)

    variance = gene_pca$sdev^2
    variance_explained = round(variance/sum(variance), 5)

    variance_df = data.frame("Variance_explained" = variance_explained,
                         "PC" = seq(1, length(variance)))

    write.table(variance_df, variance_table, sep="\t", quote=FALSE,
                row.names=FALSE)

    p_variance = ggplot(variance_df, aes(x=PC, y=Variance_explained))+
    geom_point()+
    geom_line()+
    theme_classic()+
    ylab("Variance explained (%)")+
    theme_bw() +
    theme(axis.text.x = m_text,
        axis.title.y = m_text,
        axis.title.x = m_text,
        axis.text.y = m_text,
        aspect.ratio=1)

    ggsave(variance_outfile, width=10, height=10)

    t = theme(
    axis.text.x = s_text,
    axis.title.y = m_text,
    axis.title.x = m_text,
    axis.text.y = s_text,
    legend.title = m_text,
    legend.text = m_text,
    aspect.ratio=1)

    s_c_d = scale_colour_discrete(name="Sample")
    s_s_d = scale_shape_discrete(name="Sample")

    g_p = geom_point(size=2, aes(shape=id_2, colour=id_2))

    PCs_df = data.frame(gene_pca$x)

    PCs_df$id_1 = sapply(strsplit(rownames(PCs_df), "_"), "[", 1)
    PCs_df$id_2 = sapply(strsplit(rownames(PCs_df), "_"), "[", 2)
    PCs_df$id_2 = gsub("day0", "mES Cells", PCs_df$id_2)
    PCs_df$id_2 = gsub("day", "day ", PCs_df$id_2)
    PCs_df$id_2 = relevel(as.factor(PCs_df$id_2), "mES Cells")

    write.table(PCs_df, vector_outfile,
                sep="\t", row.names=FALSE)


    p_pca1 = ggplot(PCs_df, aes(x=PC1, y=PC2)) +
    g_p + theme_bw() + t + s_c_d + s_s_d +
    xlab(paste0('PC1 (Variance explained = ' ,
                        round(100 * variance_explained[1], 1), '%)')) +
    ylab(paste0('PC2 (Variance explained = ' ,
                        round(100 * variance_explained[2], 1), '%)'))

    ggsave(pca_outfile1, width=8, height=7)

    p_pca2 = ggplot(PCs_df, aes(x=PC3, y=PC4)) +
    g_p + theme_bw() + t + s_c_d + s_s_d +
    xlab(paste0('PC3 (Variance explained = ' ,
                        round(100 * variance_explained[3], 1), '%)')) +
    ylab(paste0('PC4 (Variance explained = ' ,
                        round(100 * variance_explained[4], 1), '%)'))

    ggsave(pca_outfile2, width=8, height=7)

    loadings = data.frame(gene_pca$rotation)
    loadings = loadings[,1:10]

    write.table(loadings, pca_loadings, sep="\t", quote=FALSE)

    }''')

    df = pd.io.parsers.read_csv(infile, sep="\t", index_col=0)
    counts = Counts.Counts(df)

    counts.normalise("total-count")
    counts.log()
    # counts.transform(method="vst", design=design, blind=False)
    df = counts.table
    r_df = com.convert_to_r_dataframe(df)

    variance_outfile = P.snip(pca_loadings, "_loadings.tsv") + "_variance.png"
    variance_table = P.snip(pca_loadings, "_loadings.tsv") + "_variance.tsv"
    vector_outfile = P.snip(pca_loadings, "_loadings.tsv") + "_vectors.tsv"
    pca_outfile1 = P.snip(pca_loadings, "_loadings.tsv") + "_pca1_pca2.png"
    pca_outfile2 = P.snip(pca_loadings, "_loadings.tsv") + "_pca3_pca4.png"

    plot_PCA(r_df, variance_table, variance_outfile, vector_outfile,
             pca_outfile1, pca_outfile2, pca_loadings)


@cluster_runnable
def plotVarianceGSE65525(infiles, outfile, PCs=10):
    ''' plot varaince explained in first x PCs '''
    df = pd.DataFrame()

    for infile in infiles:
        method = P.snip(os.path.basename(infile), "_loadings.tsv").replace(
            "plots_PCA_", "").replace("dedup_", "")

        tmp_df = pd.read_table(infile.replace("loadings", "variance"), sep="\t")
        tmp_df['method'] = method
        tmp_df['cumsum'] = np.cumsum(tmp_df['Variance_explained'])

        df = df.append(tmp_df)

    df.reset_index(inplace=True)
    df.to_csv(P.snip(outfile, ".png") + ".tsv",  sep="\t", index=False)

    plotVariance = R('''
    function(df){
    library(ggplot2)

    m_text = element_text(size=20)

    t = theme(
    axis.text.x = m_text,
    axis.title.y = m_text,
    axis.title.x = m_text,
    axis.text.y = m_text,
    legend.title = m_text,
    legend.text = m_text,
    panel.grid.minor.x =element_blank(),
    panel.grid.minor.y =element_blank(),
    aspect.ratio=1)

    df$method = factor(df$method,
                       levels=c("transcriptome", "unique", "percentile",
                                "cluster", "adjacency", "directional"))

    p = ggplot(df, aes(PC, 100*cumsum, col=as.factor(method), group=method)) +
    geom_line() + theme_bw() + t +
    ylab("%% Variance Explained (Cumulative)") +
    xlab("PC") +
    scale_colour_discrete(name="Method",
                          labels=c("None", "Unique", "Percentile",
                                   "Cluster", "Adjacency", "Directional")) +
    scale_x_continuous(limits=c(1,6), breaks=seq(1,10,1)) +
    scale_y_continuous(limits=c(0,11), breaks=seq(0,11,1))

    ggsave("%(outfile)s", width=12, height=10)
    }''' % locals())

    r_df = com.convert_to_r_dataframe(df)
    plotVariance(r_df)


@cluster_runnable
def plotLoadingsGSE65525(infiles, outfile, PCs=10):
    ''' plot loadings between PCAs using different methods'''

    df = pd.DataFrame()

    for infile in infiles:
        method = P.snip(os.path.basename(infile), "_loadings.tsv").replace(
            "plots_PCA_", "").replace("dedup_", "")
        tmp_df = pd.read_table(infile, sep="\t")
        tmp_df.columns = [x + "_" + method for x in tmp_df.columns]

        df = pd.concat([df, tmp_df], axis=1)

    df.to_csv(outfile, sep="\t")


@cluster_runnable
def plotVarianceGSE53638(infile, outfile, PCs=10):
    ''' plot varaince explained in first x PCs '''
    df = pd.DataFrame()

    infiles = glob.glob(os.path.join(os.path.dirname(infile), "*_variance.tsv"))

    for infile in infiles:

        method = P.snip(os.path.basename(infile), "_variance.tsv")

        tmp_df = pd.read_table(infile.replace("loadings", "variance"), sep="\t")
        tmp_df['method'] = method
        tmp_df['cumsum'] = np.cumsum(tmp_df['Variance_explained'])

        df = df.append(tmp_df)

    df.reset_index(inplace=True)
    df.to_csv(P.snip(outfile, ".png") + ".tsv",  sep="\t", index=False)

    plotVariance = R('''
    function(df){
    library(ggplot2)

    m_text = element_text(size=20)

    t = theme(
    axis.text.x = m_text,
    axis.title.y = m_text,
    axis.title.x = m_text,
    axis.text.y = m_text,
    legend.title = m_text,
    legend.text = m_text,
    panel.grid.minor.x =element_blank(),
    panel.grid.minor.y =element_blank(),
    aspect.ratio=1)

    df$method = factor(df$method,
                       levels=c("unique", "percentile",
                                "cluster", "adjacency",
                                "directional"))

    p = ggplot(df, aes(PC, 100*cumsum, col=method, group=method)) +
    geom_line() + theme_bw() + t +
    ylab("%% Variance Explained (Cumulative)") +
    xlab("PC") +
    scale_colour_discrete(name="Method",
                          labels=c("Unique", "Percentile",
                                   "Cluster", "Adjacency", "Directional")) +
    scale_x_continuous(limits=c(1,6), breaks=seq(1,6,1)) +
    scale_y_continuous(limits=c(0,12), breaks=seq(0,12,1))

    ggsave("%(outfile)s", width=12, height=10)
    }''' % locals())

    r_df = com.convert_to_r_dataframe(df)
    plotVariance(r_df)


@cluster_runnable
def plotVectorsVsExpressionGSE53638(infile, outfile, PCs=4):
    ''' Correlate the PC vectors against the total expression '''

    print os.path.join(os.path.dirname(infile), "*_variance.tsv")
    for infile in glob.glob(os.path.join(os.path.dirname(infile), "*_variance.tsv")):

        method = P.snip(os.path.basename(infile), "_variance.tsv")

        df = pd.read_table(infile.replace("variance", "eigenvectors"), sep="\t")
        df.to_csv(outfile, sep="\t")

        r_df = com.convert_to_r_dataframe(df)

        for ix in range(0, PCs):

            PC = "PC" + str(ix + 1)

            tmp_outfile = P.snip(outfile, ".png") + "_%s_%s.png" % (PC, method)

            plotVectorExpression = R('''
            function(df){
            library(ggplot2)

            m_text = element_text(size=15)
            s_text = element_text(size=10)

            t = theme(
            axis.text.x = s_text,
            axis.title.y = m_text,
            axis.title.x = m_text,
            axis.text.y = s_text,
            legend.title = m_text,
            legend.text = m_text,
            panel.grid.minor.x =element_blank(),
            panel.grid.minor.y =element_blank(),
            aspect.ratio=1)

            p = ggplot(df, aes(id_expression, %(PC)s, col=id_2, group=id_2)) +
            geom_point() + t +
            scale_colour_discrete(name="Stage") +
            facet_grid(~id_2) +
            geom_smooth(method="lm")


            ggsave("%(tmp_outfile)s", width=12, height=5)

            }''' % locals())

            print plotVectorExpression

            plotVectorExpression(df)


@cluster_runnable
def plotFacettedEditPlots(infiles, outfile):
    ''' combine the edit distance from all samples and plot the edit
    distances '''

    def sample2experiment(sample):
        
        sample2day ={"SRR1058003": "Day 0",
                     "SRR1058023": "Day 14",
                     "SRR1058032": None,
                     "SRR1058038": None,
                     "SRR1784310": "mES Cells",
                     "SRR1784313": "Day 2",
                     "SRR1784314": "Day 4",
                     "SRR1784315": "Day 7"}

        if sample in ["SRR1058003", "SRR1058023",
                      "SRR1058032", "SRR1058038"]:
            experiment = "SCRB-Seq"

        elif sample in ["SRR1784310", "SRR1784313",
                        "SRR1784314", "SRR1784315"]:
            experiment = "inDrop-Seq"

        day = sample2day[sample]

        return day, experiment


    df = pd.DataFrame()
    keep_columns = ["post", "post_null", "edit_distance", "method"]
    for infile in infiles:
        sample = os.path.basename(infile).replace(
            "_merged_edit_distances.tsv", "")
        tmp_df = pd.read_table(infile, usecols=keep_columns)
        day, experiment = sample2experiment(sample)
        if day:
            tmp_df['sample'] = day
            tmp_df['experiment'] = experiment
            df = pd.concat((df, tmp_df))

    df.to_csv(outfile, sep="\t")
    df = df[df['method'].isin(["Unique", "Percentile", "Directional"])]

    plotFacettedEditDistance = R('''
    function(df, samples, plotfile){
    library(ggplot2)
    library(reshape)
    library(plyr)
    library(grid)

    df = df[df$edit_distance!="Single_UMI",]
    df$edit_distance = as.numeric(as.character(df$edit_distance))
    df_null = df[df$method=="Unique",]
    df_null$method="Null"
    df_null$post = df_null$post_null
    df = rbind(df, df_null)

    max_distance = max(df$edit_distance)
    col_range = topo.colors(max_distance, alpha = 0.9)

    df = df[df$edit_distance<max_distance,]
    df = df[order(df$edit_distance),]

    l_txt = element_text(size=35)

    df$method = factor(df$method,
                       levels=c("Unique", "Percentile", "Cluster",
                                "Adjacency", "Directional",
                                "Null"))

    p = ggplot(df, aes(y=post, x=method)) +

        geom_bar(aes(fill=as.factor(edit_distance), colour="grey50"),
                     stat="identity", position="fill", width=1) +

        xlab("") +

        ylab("Fraction") +

        scale_fill_manual(name="Edit distance", values=col_range) +

        scale_colour_manual(guide="none", values="grey60") +

        theme_bw() +

        theme(axis.text.y=l_txt,
              axis.text.x=element_text(size=35, angle=90,
                                       hjust=1, vjust=0.5),
              axis.title.x=l_txt, axis.title.y=l_txt,
              strip.background = element_rect(colour="grey60", fill="grey90"),
              legend.text=l_txt, strip.text=l_txt,
              legend.title=l_txt,
              panel.margin = unit(1, "lines"),
              legend.position="right",
              aspect.ratio=1.5,
              legend.key.size=unit(2, "line"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +

        facet_grid(experiment~sample) +

    ggsave(plotfile, width=((5*samples) + 2), height=7.5)
    }''')

    for experiment in set(df['experiment'].tolist()):
        tmp_df = df[df['experiment'] == experiment]
        samples = set(tmp_df['sample'].tolist())
        plotfile = P.snip(outfile, ".tsv") + "_%s.png" % experiment
        plotFacettedEditDistance(tmp_df, len(samples), plotfile)
