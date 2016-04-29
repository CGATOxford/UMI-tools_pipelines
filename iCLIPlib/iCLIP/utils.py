import re
import numpy as np
import pandas as pd
import CGAT.GTF as GTF


AMBIGUITY_CODES = {'M': 'AC',
                   'R': 'AG',
                   'W': 'AT',
                   'S': 'CG',
                   'Y': 'CT',
                   'K': 'GT',
                   'V': 'ACG',
                   'H': 'ACT',
                   'D': 'AGT',
                   'B': 'CGT',
                   'N': 'CGAT'}


##################################################
def IUPAC2Regex(sequence):

    for code, regex in AMBIGUITY_CODES.iteritems():
        sequence = re.sub(code, '[%s]' % regex, sequence)

    return sequence


##################################################
class TranscriptCoordInterconverter:
    ''' A class to interconvert between genome co-ordinates
    and transcript co-ordinates. Implemented as a class because
    there are expected to be many calls against the same transcript,
    so time can be saved by precomputation

    TranscriptCoordInterconverter.genome2transcript should be the
    interverse of TranscriptCoordInterconverter.transcript2genome.

    That is

    if myConverter = TranscriptCoordInterverter(transcript)
    
    then
    
    myConverter.genome2transcript(myConverter.transcript2genome(x)) == x
    
    and

    myConverter.transcript2genome(myConverter.genome2transcript(x)) == x'''

    def __init__(self, transcript, introns=False):
        ''' Pre compute the conversions for each exon '''

        if not introns:
            intervals = GTF.asRanges(transcript, feature="exon")
        else:
            intervals = GTF.toIntronIntervals(transcript)
        
        # get strand
        self.strand = transcript[0].strand

        # store transcript_id
        self.transcript_id = transcript[0].transcript_id

        # sort the exons into "transcript" order
        if self.strand == "-":
            intervals.sort(reverse=True)
            intervals = [(y-1, x-1) for x, y in intervals]
        else:
            intervals.sort(reverse=False)

        self.offset = intervals[0][0]
        self.genome_intervals = [map(abs, (x-self.offset, y-self.offset))
                                 for x, y in intervals]

        interval_sizes = [abs(y-x) for x, y in intervals]

        total = 0
        transcript_intervals = [None]*len(interval_sizes)

        for i in range(len(interval_sizes)):
            transcript_intervals[i] = (total,
                                       interval_sizes[i] + total)
            total += interval_sizes[i]
        
        self.transcript_intervals = transcript_intervals
        self.length = transcript_intervals[-1][1]

    def genome2transcript(self, pos):
        ''' Convert genome coordinate into transcript coordinates.
        pos can be a single value or a nunpy array like object.
        Passing an array ensures that the transcript is only
        searched once, ensuring O(n) performance rather than
        O(nlogn)'''

        if len(pos) == 0:
            return np.array([])

        try:
            relative_pos = pos - self.offset
        except TypeError:
            relative_pos = np.array(pos) - self.offset
        
        if self.strand == "-":
            relative_pos = relative_pos * -1

        ordering = np.argsort(relative_pos)
        relative_pos = np.sort(relative_pos)

        # pre allocate results list for speed
        try:
            results = np.zeros(len(relative_pos))
        except TypeError:
            relative_pos = np.array([relative_pos])
            results = np.zeros(1)

        i = 0
        i_max = len(relative_pos)

        # do linear search for correct exon
        for exon, interval in enumerate(self.genome_intervals):

            if relative_pos[i] < interval[0]:

                raise ValueError("Position %i is not in transcript %s" %
                                 (pos[i], self.transcript_id))
            
            while relative_pos[i] < interval[1]:
                
                pos_within_exon = relative_pos[i]-interval[0]
                transcript_exon = self.transcript_intervals[exon]
                transcript_position = transcript_exon[0] + pos_within_exon
                
                results[i] = transcript_position
                i += 1
                if i == i_max:
                    return results[ordering]

        # exon has not been found
       
        raise ValueError("Position %i (%i relative) is not in transcript %s\n exons are %s" %
                         (pos[i], relative_pos[i], self.transcript_id, self.genome_intervals))

    def transcript2genome(self, pos):
        ''' Convert transcript coodinate into genome coordinate,
        pos can be a single value or a nunpy array like object.
        Passing an array ensures that the transcript is only
        searched once, ensuring O(n) performance rather than
        O(nlogn)'''
    
        try:
            if len(pos) == 0:
                return np.array([])
        except TypeError:
            pos = np.array([pos])

        # Converting a list is only efficient if the list is ordered
        # however want to be able to return list in the same order it
        # arrived, so remember the order and then sort.
        ordering = np.argsort(pos)
        pos = np.sort(pos)

        # pre allocate results list for speed
       
        results = np.zeros(len(pos))
        
        i = 0
        i_max = len(pos)

        # do linear search for correct exon
        for exon, interval in enumerate(self.transcript_intervals):

            while pos[i] < interval[1]:
                pos_within_exon = pos[i] - interval[0]
                genome_exon = self.genome_intervals[exon]
                relative_genome_position = genome_exon[0] + pos_within_exon

                if self.strand == "-":
                    results[i] = (self.offset - relative_genome_position)
                    i += 1
                else:
                    results[i] = (self.offset + relative_genome_position)
                    i += 1

                if i == i_max:
                    return results[ordering]
  
        # beyond the end of the transcript
        ValueError("Transcript postion %i outside of transcript %s" %
                   (pos[i], self.transcript_id))

    def transcript_interval2genome_intervals(self, interval):
        '''Take an interval in transcript coordinates and returns
        a list of intervals in genome coordinates representing the
        interval on the genome '''

        outlist = []
        for exon in self.transcript_intervals:
            if interval[0] < exon[1]:
                start = interval[0]
                if interval[1] <= exon[1]:
                    outlist.append((start, interval[1]))
                    break
                else:
                    outlist.append((start, exon[1]))
                    interval = (exon[1], interval[1])
       
        genome_list = [tuple(self.transcript2genome((x, y-1))) for
                       x, y in outlist]
        
        # these intervals are zero based-closed. Need to make half open

        if self.strand == "+":
            genome_list = [(x, y+1) for x, y in genome_list]
        else:
            genome_list = [(y, x+1) for x, y in genome_list]

        return sorted(genome_list)


##################################################
def randomiseSites(profile, start, end, keep_dist=True):
    '''Randomise clipped sites within an interval (between start and end)
    if keep_dist is true, then reads on the same base are kept togehter'''

    if keep_dist:

        profile = profile.copy()
        profile.index = np.random.choice(
            np.arange(start, end), profile.size, replace=False)
        profile = profile.sort_index()
        return profile

    else:
        randomised = np.random.choice(
            np.arange(start, end), profile.sum(), replace=True)
        idx, counts = np.unique(randomised, return_counts=True)
        randomised = pd.Series(counts, index=idx).sort_index()
        return randomised


##################################################
def spread(profile, bases, reindex=True, right_bases=None):
       
    if right_bases:
        window = bases+right_bases
    else:
        window = 2 * bases + 1
        right_bases = bases

    if reindex:
        start = int(profile.index.min() - window)
        end = int(profile.index.max() + window+1)
        profile = profile.reindex(np.arange(start, end), fill_value=0)

    result = pd.rolling_sum(profile, window=window, center=False).dropna()
    result.index = result.index - (right_bases + 1)
    return result


##################################################
def rand_apply(profile, exon, n, func, keep_dist=False,
               *args, **kwargs):
    '''Randomise a profile multiple times and apply a function
    to each randomised profile.

        :param profile: a profile with the number of reads at each base
        :type profile: pandas.Series
        :param exon: a GTF entry specifiying the boundaries to randomiseSites
                     between.
        :type exon: CGAT.GTF.Entry
        :param n: The number of randomisations to apply
        :param func: The function to apply to each randomisation
        :param keep_dist: Keep the read-per-base distribution constant

        :rtype: pandas.Series or pandas.DataFrame

    '''

    dummy = pd.Series(range(n))
 
    def _inner_func(x):
        rand = randomiseSites(profile, exon.start, exon.end,
                              keep_dist=keep_dist)
        return func(rand, *args, **kwargs)

    return dummy.apply(_inner_func)


