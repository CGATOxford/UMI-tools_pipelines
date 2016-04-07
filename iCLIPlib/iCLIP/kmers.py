import re
import pandas as pd
import numpy as np
import collections
import itertools

import CGAT.GTF as GTF
import CGAT.Experiment as E

from utils import spread, rand_apply, AMBIGUITY_CODES
from counting import count_transcript, count_intervals


##################################################
def find_all_matches(sequence, regexes):

    def _find_regex(regex):
        matches = regex.finditer(sequence)
        mStarts = [m.start() for m in matches]
        return np.asarray(mStarts)

    return regexes.apply(_find_regex)


##################################################
def pentamer_frequency(profile, length, regex_matches, nSpread=15):
    '''Calculate the frequency of the each of a collection
    of sequence regexes on the provided read profile, and the
    coresponding sequence

        :param profile: A profile of the number of reads at each base
        :type profile: pandas.Series
        :param length: Length of the sequence represented by profile
        :param regex_matches: A pandas.Series of the locations of hits
                              for a set of regexs, as returned by 
                              find_all_matches
        :param nSpread: How far either side of each read to consider

        :rtype: pandas.Series with the count of each regex'''

    kmer = len(regex_matches.index.values[0])
    profile = profile.reindex(
        np.arange(-nSpread, length+nSpread-kmer), fill_value=0)
    profile = spread(profile, nSpread, False, nSpread - kmer)
    profile = profile.values

    def _count_regex(hits):

        if hits.size:
            return profile[hits].sum()
        else:
            return 0

    return regex_matches.map(_count_regex)
    

##################################################
LiteExon = collections.namedtuple('LiteExon', "start, end")


##################################################
def pentamer_enrichment(gtf_chunk_iterator, bam, fasta, kmer_length=5,
                        randomisations=100,
                        seperate_UTRs=False, spread=15,
                        pool=None):
    ''' This function calculates the enrichment of pentamers around
    CLIP'd sites

        :type gtf_chunk_iterator: an iterator that returns lists of CGAT.GTF
                                  Entries.
        :type bam: pysam.AlignmentFile
        :param seperate_UTRs: Treat UTRs as seperate areas from the rest of the
                              gene. Requires GTF file to contain CDS entries.
        :param spread: Number of bp around each base to use for sequence.
        :type fasta: CGAT.IndexedFasta
        :param pool: If present work will be parallelize across the worker pool
        :type pool: multiprocessing.Pool

        :rtype: pandas.Series with z-values for each pentamer
    '''

    bases = AMBIGUITY_CODES.keys()
    bases.remove("N")
    kmers = itertools.product('CGAT', repeat=kmer_length)
    kmers = ["".join(kmer) for kmer in kmers]
    regexs = [re.compile(kmer) for kmer in kmers]
    regexs = pd.Series(regexs, index=kmers)
    
    observed_kmer_counts = pd.Series(0, index=kmers)
    randomised_kmer_counts = pd.DataFrame(0,
                                          index=np.arange(randomisations),
                                          columns=kmers)

    args = ((profile, sequence,
             regexs, spread, randomisations)
            for profile, sequence in
            _get_counts_and_sequence(gtf_chunk_iterator, bam, fasta))

    if pool:
        results_iterator = pool.imap_unordered(_par_call, args)
    else:
        results_iterator = (_get_regex_frequencies(*arg) for arg in args)

    for observed, rands in results_iterator:

        observed_kmer_counts += observed
        randomised_kmer_counts += rands

    randomised_kmer_counts = randomised_kmer_counts.transpose()
    means = randomised_kmer_counts.mean(axis=1)
    sd = randomised_kmer_counts.std(axis=1)

    means.name = "mean"
    sd.name = "sd"
    observed_kmer_counts.name = "count"

    results = pd.concat([observed_kmer_counts, means, sd], axis=1)
    return results.apply(lambda x: (x['count'] - x['mean'])/x['sd'], 1)


def _get_counts_and_sequence(gtf_iterator, bam, fasta,
                             seperate_UTRs=False):
    '''Called by pentamer_enrichment. This function will return an iterator
    that yeilds tuples of profiles accross transcripts or introns and the
    sequence for which the profile is determined'''

    for transcript in gtf_iterator:

        E.debug("Counting transcript %s" % transcript[0].transcript_id)
        contig, strand = transcript[0].contig, transcript[0].strand

        # exons
        exons = GTF.asRanges(transcript, "exon")
        sequence = "".join(fasta.getSequence(contig, strand, exon[0], exon[1])
                           for exon in exons)
        exon_counts = count_transcript(transcript, bam)
        yield (exon_counts, sequence)

        # introns
        intron_intervals = GTF.toIntronIntervals(transcript)
        intron_counts = count_intervals(bam, intron_intervals, contig, strand)

        if intron_counts.sum() == 0:
            continue

        for intron in intron_intervals:
            
            seq = fasta.getSequence(contig, strand, intron[0], intron[1])
            profile = intron_counts.loc[float(intron[0]):float(intron[1])]
            profile.index = profile.index - intron[0]
            yield (profile, seq)

                
def _get_regex_frequencies(profile, sequence, regexs,
                           spread, randomisations):
    '''Called by pentamer_enrichment to get the frequencies
    counts for observed and randomisations. Allows
    parallelisation'''

    regex_matches = find_all_matches(sequence, regexs)
    length = len(sequence)
    boundaries = LiteExon(0, length)
    observed_kmer_counts = pentamer_frequency(profile,
                                              length,
                                              regex_matches,
                                              spread)
    randomised_kmer_counts = rand_apply(profile, boundaries,
                                        randomisations,
                                        pentamer_frequency,
                                        length=length,
                                        regex_matches=regex_matches,
                                        nSpread=spread)

    return observed_kmer_counts, randomised_kmer_counts


def _par_call(args):
    return _get_regex_frequencies(*args)
