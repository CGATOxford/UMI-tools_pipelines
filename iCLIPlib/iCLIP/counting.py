''' These functions deal with calculating the read counts within a genomic region. 
in general they return a pandas Series of read counts. '''


import pandas as pd
import collections

import CGAT.Experiment as E
import CGAT.GTF as GTF

from utils import TranscriptCoordInterconverter


def find_first_deletion(cigar):
    '''Find the position of the the first deletion in a
    read from the cigar string, will return 0 if no deletion
    found '''

    position = 0
    for operation, length in cigar:

        if operation == 2:
            return position
        else:
            position += length

    return position


##################################################
def getCrosslink(read):
    ''' Finds the crosslinked base from a pysam read.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.

    returns a tuple with the position of the read and one of the following
    categories:

        * truncated_neg
  
        * truncated_pos
      
        * deletion_neg

        * deletion_pos


    to record whether the position came from a truncation or a deletion '''

    if 'D' not in read.cigarstring:
        if read.is_reverse:
            pos = read.aend

        else:
            pos = read.pos - 1

    else:
        if read.is_reverse:
            cigar = reversed(read.cigar)
            position = find_first_deletion(cigar)
            pos = read.aend - position - 1

        else:
            position = find_first_deletion(read.cigar)
            pos = read.pos + position

    return pos


##################################################
def countChr(reads, chr_len, dtype='uint16'):
    ''' Counts the crosslinked bases for each read in the pysam rowiterator
    reads and saves them in pandas Series: those on the positive strand
    and those on the negative strand. The Series are indexed on genome position,
    and are sparse.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.
    
    The dtype to use internally for storage can be specified. Large types
    reduce the chance of overflow, but require more memory. With 'uint16'
    the largest count that can be handled is 255. Data is stored sparse,
    so memory is less of a problem. Overflow will cause a ValueError.

    returns a tuple of pandas Series objects, with the positive and negative
    strand arrays and also a counter object that contains the counts for each
    type of site. '''

    pos_depths = collections.defaultdict(int)
    neg_depths = collections.defaultdict(int)

    counter = 0

    for read in reads:
        
        pos = getCrosslink(read)
        counter += 1

        if read.is_reverse:
            neg_depths[float(pos)] += 1
        else:
            pos_depths[float(pos)] += 1

    try:
        pos_depths = pd.Series(pos_depths, dtype=dtype)
    except ValueError:
        pos_depths = pd.Series({}, dtype=dtype)

    try:
        neg_depths = pd.Series(neg_depths, dtype=dtype)
    except ValueError:
        neg_depths = pd.Series({}, dtype=dtype)
 
    # check for integer overflow: counter sum should add up to array sum
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter == array_sum:
        raise (ValueError,
               "Sum of depths is not equal to number of "
               "reads counted, possibly dtype %s not large enough" % dtype)
    
#    E.debug("Counted %i truncated on positive strand, %i on negative"
#            % (counter.truncated_pos, counter.truncated_neg))
#    E.debug("and %i deletion reads on positive strand, %i on negative"
#            % (counter.deletion_pos, counter.deletion_neg))
    
    return (pos_depths, neg_depths, counter)


##################################################
def count_intervals(bam, intervals, contig, strand=".", dtype='uint16'):
    ''' Count the crosslinked bases accross a transcript '''

    chr_len = bam.lengths[bam.gettid(contig)]
    exon_counts = []
    for exon in intervals:
        
        # X-linked position is first base before read: need to pull back
        # reads that might be one base out. Extra bases will be filtered out
        # later.
        try:
            reads = bam.fetch(reference=contig,
                              start=max(0, exon[0]-1),
                              end=exon[1]+1)
        except ValueError as e:
            E.debug(e)
            E.warning("Skipping intervals on contig %s as not present in bam"
                      % contig)
            return pd.Series()

        count_results = countChr(reads, chr_len, dtype)

        # fetch pulls back any reads that *overlap* the specified coordinates
        # exlude Xlinked bases outside the interval (prevents double counting)

        if strand == "+":
            if len(count_results[0]) > 0:
                exon_counts.append(count_results[0].sort_index().loc[
                    float(exon[0]):float(exon[1]-1)])
        elif strand == "-":
            if len(count_results[1]) > 0:
                exon_counts.append(count_results[1].sort_index().loc[
                    float(exon[0]):float(exon[1]-1)])
            
        else:
            sum_counts = count_results[0].loc[
                (count_results[0].index >= exon[0]) &
                (count_results[0].index < exon[1])].add(
                    count_results[1].loc[(count_results[1].index >= exon[0]) &
                                         (count_results[1].index < exon[1])],
                    fill_value=0)
            exon_counts.append(sum_counts)

    if len(exon_counts) == 0:
        transcript_counts = pd.Series()
    else:
        transcript_counts = pd.concat(exon_counts)
    # transcript_counts = transcript_counts.sort_index()

    return transcript_counts


##################################################
def count_transcript(transcript, bam, flanks=0):
    '''Count clip sites from a bam and return a Series which transcript 
    domain coordinates.

        :type transcript: list of CGAT.GTF.Entry or similar
        :type bam: pysam.AlignmentFile
        :param flanks: If specified, flanks of this length are counted
                       and returned to the flank3 and flank5 section
                       of the now multiindexed return series
        :rtype: pandas.Series

    '''

    import pandas
    exons = GTF.asRanges(transcript, "exon")
    
    counts = count_intervals(bam, exons,
                             contig=transcript[0].contig,
                             strand=transcript[0].strand)
    
    coords_translator = TranscriptCoordInterconverter(transcript)
    
    counts.index = coords_translator.genome2transcript(counts.index)

    if flanks == 0:
        return counts
    else:

        if counts.sum() > 0:
            counts.index = pandas.MultiIndex.from_tuples(
                zip(*(['exons']*len(counts), counts.index)),
                names=["region", "base"])
        else:
            # deal with empty case
            counts.index = pandas.MultiIndex(
                levels=[["flank5", "exon", "flank3"],
                        []],
                labels=[[], []],
                names=["region", "base"])
 
        transcript_min = min([a for a, b in exons])
        transcript_max = max([b for a, b in exons])

        E.debug("Counting 5 flank: (%i, %i)" % (transcript_min - flanks,
                                                transcript_min))
        flank5_counts = count_intervals(bam, [(transcript_min - flanks,
                                              transcript_min)],
                                        contig=transcript[0].contig,
                                        strand=transcript[0].strand)
        E.debug(flank5_counts)
        E.debug("Counting 3 flank: (%i, %i)" % (transcript_max,
                                                transcript_max + flanks))
        flank3_counts = count_intervals(bam, [(transcript_max,
                                              transcript_max + flanks)],
                                        contig=transcript[0].contig,
                                        strand=transcript[0].strand)
        E.debug(flank3_counts)

        if transcript[0].strand == "-":
            index5 = ['flank3']
            index3 = ['flank5']
            bases5 = transcript_min - flank5_counts.index.values
            bases3 = transcript_max + flanks - flank3_counts.index.values

        else:
            index5 = ['flank5']
            index3 = ['flank3']
            bases5 = flank5_counts.index.values - transcript_min + flanks
            bases3 = flank3_counts.index.values - transcript_max

        if flank3_counts.sum() > 0:
            flank3_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index3*len(bases3), bases3)), names=["region", "base"])

        if flank5_counts.sum() > 0:
            flank5_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index5*len(bases5), bases5)), names=["region", "base"])

        E.debug("After direction detection and reindexing:")
        E.debug(flank3_counts)
        E.debug(flank5_counts)

        result = pandas.concat([flank5_counts, counts, flank3_counts])
        return result

