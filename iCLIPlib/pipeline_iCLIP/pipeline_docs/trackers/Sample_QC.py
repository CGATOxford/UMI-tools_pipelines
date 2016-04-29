
from iCLIPTracker import *
import numpy as np
import scipy as sp
from collections import OrderedDict

class ContextStats(iCLIPTracker):

    categories = {"antisense": "Long ncRNA",
                  "intron": "Intron",
                  "retained_intron": "Retained intron",
                  "lincRNA": "Long ncRNA",
                  "miRNA": "Small ncRNA",
                  "miscRNA": "Other",
                  "none": "Genomic",
                  "nonsense_mediated_decay": "Other",
                  "polymorphic_pseudogene": "Long ncRNA",
                  "processed_pseudogene": "Long ncRNA",
                  "processed_transcript": "Other",
                  "protein_coding": "Protein coding exon",
                  "pseudogene": "Long ncRNA",
                  "rRNA": "Ribosomal RNA",
                  "snoRNA": "Small ncRNA",
                  "snRNA": "Small ncRNA",
                  "unitary_pseudogene": "Long ncRNA",
                  "unprocessed_pseudogene": "Long ncRNA"}

#   slices = ["mapping", "deduped", "clusters"]
    method = "deduped"
    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM %(method)s_context_stats")

    def __call__(self, track, slice=None):

        statement = ''' SELECT category, alignments
                      FROM %(method)s_context_stats
                      WHERE track = '%(track)s' AND category != 'total' '''

        results = self.getDataFrame(statement)
        results.category = [self.categories[x] if x in self.categories else "Other"
                            for x in results.category]
        results = results.groupby("category").sum()
                            
        results = results.groupby(level="category").sum().reset_index()

        return results


class ContextRepresentation(ContextStats):


    categories = ["antisense",
                  "intron",
                  "lincRNA",
                  "miRNA",
                  "miscRNA",
                  "none",
                  "nonsense_mediated_decay",
                  "polymorphic_pseudogene",
                  "processed_pseudogene",
                  "processed_transcript",
                  "protein_coding",
                  "pseudogene",
                  "rRNA",
                  "snoRNA",
                  "snRNA",
                  "unitary_pseudogene",
                  "unprocessed_pseudogene"]

    def __call__(self, track):

        slice = "deduped"
        categories = "','".join(self.categories)
        statement = ''' SELECT category, 
                               alignments,
                               nbases,
                               (alignments + 0.0) / (SELECT SUM(alignments)
                                  FROM %(slice)s_context_stats
                                  WHERE track = '%(track)s'
                                       AND category != 'total') as precent_alignments,
                               ( nbases + 0.0) /
                               (SELECT SUM(nbases) FROM reference_context_interval_stats)
                                 AS percent_bases
                      FROM %(slice)s_context_stats cstats
                       INNER JOIN
                            reference_context_interval_stats as ints
                          ON ints.track = cstats.category
                      WHERE cstats.track = '%(track)s'
                         AND cstats.category IN ('%(categories)s')  
                      '''

       
        results = self.getAll(statement)

        return results

        
class ContextSaturation(ProjectTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT category FROM saturation_context_stats")

    def getSlices(self):
        return self.getValues("SELECT DISTINCT track FROM saturation_context_stats")

    def __call__(self, track, slice):

        factor, replicate = re.match("(.+)\-.+\-(.+)", slice).groups()

        statement = '''SELECT subset, alignments, 
                       '%(factor)s' as factor, '%(replicate)s' as replicate
                      FROM saturation_context_stats
                      WHERE track = '%(slice)s' AND
                            category = '%(track)s' '''

        results = self.getAll(statement)
        total = max(results["alignments"])
        results["alignments"] = [float(x)/total for x in results["alignments"]]
        return results


class AlignmentSaturation(ProjectTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM saturation_context_stats")

    def __call__(self, track):

        factor, replicate = re.match("(.+)\-.+\-(.+)", track).groups()

        statement = '''SELECT subset, counts, 
                       '%(factor)s' as factor, '%(replicate)s' as replicate
                      FROM subset_bam_stats
                      WHERE track = '%(track)s' AND
                            category = 'reads_total' '''

        results = self.getAll(statement)
        return results

def saturation_function(s, n):
    return n*(1.0-(1.0-(1.0/n))**s)

def saturation_function2(s, vmax, km):
    return s*vmax/(km + s)


class LibrarySizes(ProjectTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM saturation_context_stats")

    def __call__(self, track):
        
        statement = '''SELECT subset, counts as alignments
                      FROM subset_bam_stats
                      WHERE track = '%(track)s' AND
                      category = 'alignments_total' '''
        
        mapper = PARAMS["mappers"]
        results = self.getAll(statement)
        total = self.getValue(''' SELECT sum(reads_total) FROM mapping.view_mapping
                                  WHERE track LIKE '%(track)s_%%.%(mapper)s' ''')


        xdata = np.array(results["subset"])
        xdata = xdata*total
        ydata = np.array(results["alignments"])
    
        print xdata
        print ydata


        popt = sp.optimize.curve_fit(globals()[self.function],
                                     xdata, ydata, (total,)+self.p0)
        popt = popt[0]
        total_unique_sequenced = max(ydata)
        percent_saturation = total_unique_sequenced/popt[0]

        print popt
        results["expected_unique"] = [globals()[self.function](s, *popt)
                                      for s in xdata]
        results["library_size"] = [popt[0] for s in xdata]
        results["subset"] = xdata
        #P(X>0) = 1 - P(X=0) = 1- (n 0) * p**0 * q**n = 1 - q ** n
        # n required for 95% capture:
        q = 1 - 1/popt[0]
        n_required = np.log(0.05)/np.log(q)

        return results

class LibrarySize_Binom(LibrarySizes):

    function = "saturation_function"
    p0=tuple()

class LibrarySize_mm(LibrarySizes):
    function = "saturation_function2"
    p0=(1,)

class fit_stats(LibrarySizes):

    def __call__(self, track):

        results = LibrarySizes.__call__(self, track)
        lib_size = results["library_size"][0]
        total_unique_sequenced = max(results["alignments"])
        percent_saturation = total_unique_sequenced/lib_size
        
        #P(X>0) = 1 - P(X=0) = 1- (n 0) * p**0 * q**n = 1 - q ** n
        # n required for 80% capture:
        q = 1 - 1/lib_size
        n_required = np.log(0.05)/np.log(q)

        mapper = PARAMS["mappers"]
        sequenced = self.getValue(''' SELECT sum(reads_total) FROM mapping.view_mapping
                                  WHERE track LIKE '%(track)s_%%.%(mapper)s' ''')

        return {"Library Size": lib_size,
                "Percent Saturation": percent_saturation*100,
                }


class mm_fit_stats(LibrarySize_mm, fit_stats):
    pass

        
class UMI_stats(ProjectTracker):

    pattern = "(.+[^_])umi_stats"

    def __call__(self, track):

        statement = '''SELECT track as sample, UMI, (Count +0.0)/sum_count as freq
                       FROM %(track)sumi_stats as umi
                       INNER JOIN
                        sample_table as samples 
                     ON umi.Sample = samples.barcode
                       INNER JOIN
                        (SELECT Sample, sum(Count) as sum_count
                         FROM %(track)sumi_stats
                         GROUP BY Sample) as sum_stats
                      ON umi.Sample = sum_stats.Sample
                       '''

        return self.getAll(statement)


class ReadsPerSample(ProjectTracker):

    pattern = "(.+[^_])umi_stats"

    def __call__(self, track):
        
        mapper = PARAMS["mappers"]
        statement = '''SELECT samples.track as sample,
                              sum(count) as total
                       FROM %(track)sumi_stats as umi
                       INNER JOIN
                        sample_table as samples
                       ON samples.barcode = umi.Sample
          
                       GROUP BY Sample '''

        return self.getAll(statement)



class PercentDemuxed(ProjectTracker):

    pattern = "(.+[^_])umi_stats"

    def __call__(self, track):
        print "called with track %s" % track
        mapper = PARAMS["mappers"]
        statement = '''SELECT samples.track as sample,
                              (vm.reads_total +0.0)/sum(count) as demuxed
                       FROM %(track)sumi_stats as umi
                       INNER JOIN
                        sample_table as samples
                       ON samples.barcode = umi.Sample
                       INNER JOIN
                         mapping.view_mapping as vm
                       ON vm.track =  samples.track || '_%(track)s' || '.%(mapper)s'
                       GROUP BY Sample '''

        return self.getAll(statement)


class PercentMapped(ProjectTracker):

    tracks = ["merged"]
    def __call__(self, track):

        mapper = PARAMS["mappers"]
        statement = '''SELECT samples.track as sample,
                              (vm.reads_mapped +0.0)/vm.reads_total as mapped
                       FROM
                        sample_table as samples
                       INNER JOIN
                         mapping.view_mapping as vm
                       ON vm.track = "%(track)s_" || samples.track || '.%(mapper)s'
                       GROUP BY Sample '''
        return self.getAll(statement)


class PercentDeDuped(ProjectTracker):

    tracks = ["merged"]
    def __call__(self, track):

        mapper=PARAMS["mappers"]

        statement = '''SELECT dd.track as sample,
                       (dd.counts + 0.0)/vm.reads_mapped as p_unique
                       FROM
                        deduped_bam_stats as dd
                       INNER JOIN
                        mapping.view_mapping as vm
                       ON vm.track = '%(track)s_' || dd.track || '.%(mapper)s'
                       WHERE dd.category = 'reads_mapped' '''

        return self.getAll(statement)


class FinalReads(ProjectTracker):

    def __call__(self, track):

        statement = '''SELECT dd.track as sample,
                       dd.counts as reads_mapped
                       FROM
                        deduped_bam_stats as dd
                       
                       WHERE dd.category = 'reads_mapped' '''

        return self.getAll(statement)

class PercentSpliced(ProjectTracker):
    ''' Percent of deduped reads that are spliced '''

    def __call__(self, track):

        statement = ''' SELECT ns.track as Track,
                               (ns.nspliced + 0.0)/dbs.counts as pspliced
                        FROM deduped_nspliced as ns
                          INNER JOIN 
                            deduped_bam_stats as dbs
                          ON ns.track = dbs.track 
                       WHERE category = 'alignments_mapped' '''
        return self.getAll(statement)

class Reproducibility(ProjectTracker):

    table = "experiment_reproducibility"
    pattern = "(.+\-.+)\-.+\.bam"

    def getTracks(self):

        tracks = self.getValues("SELECT DISTINCT Experiment FROM %(table)s")
        return tracks

    def getSlices(self):

        return self.getValues("SELECT DISTINCT fold FROM %(table)s")

    def __call__(self, track, slice):

        statement = ''' SELECT  Track as sample, level, (hits+0.0)/totals as reproducibility
                        FROM %(table)s
                        WHERE Experiment = '%(track)s' AND fold = %(slice)s 
                         AND totals > 100'''

        results =  self.getAll(statement)
        results['Replicate'] = [re.match(".+\-.+\-(.+).bam", t).groups()[0] for t in results["sample"]]
        return results

class NormReproducibility(Reproducibility):

    def __call__(self, track, slice):
        
        statement = ''' SELECT Track, level, hits, totals
                        FROM %(table)s
                        WHERE fold = %(slice)s AND Experiment = '%(track)s'
                              AND totals > 100 '''

        results = self.getAll(statement)

        replicate = [re.match(".+\-.+\-(.+).bam", t).groups()[0] for t in results["Track"]]

        totals = self.getValues('''
           SELECT totals
           FROM %(table)s
           WHERE level = 1 AND fold = 1 AND Experiment = '%(track)s' ''')
        print totals
        totals = map(int, totals)
        print totals
        results["reproducibility"] = []

        if slice == 1:
            totals = sum(totals)
            for i in range(len(results["Track"])):
                total_of_others = totals - results["totals"][i]
                total = min(results["totals"][i], total_of_others)
                repro = float(results["hits"][i])/total
                results["reproducibility"].append(repro)
        elif slice == 2:
            for i in range(len(results["Track"])):
                total = min(list(totals) + [results["totals"][i]])
                repro = float(results["hits"][i])/total
                results["reproducibility"].append(repro)

        return {"Replicate": replicate,
                "level": results["level"],
                "reproducibility": results["reproducibility"]}
                

class ReproducibilityAll(Reproducibility):

    table = "all_reproducibility"

    def getTracks(self):

        tracks = self.getValues("SELECT DISTINCT Track FROM %(table)s")
        return tracks

    def __call__(self, track, slice):

        statement = ''' SELECT  Track as sample, level, (hits+0.0)/totals as reproducibility
                        FROM %(table)s
                        WHERE Track = '%(track)s' AND fold = %(slice)s 
                         AND totals > 100'''

        results = self.getAll(statement)
        results['Replicate'] = [re.match(".+\-.+\-(.+).bam", t).groups()[0] for t in results["sample"]]
        return results


class ReproducibilityVsControl(Reproducibility):
    table = "reproducibility_vs_control"


class ReproducibilityReplicateVsControl(Reproducibility):
    '''measure of interest is the ratio of the reproducibility in
    replicates and controls '''


    slices = ["1","All"]

    def getTracks(self):

        tracks = self.getValues("SELECT DISTINCT Experiment FROM reproducibility_vs_control")
        return tracks

    def __call__(self, track, slice):

        if slice == "1":
            rep_fold = 1
            control_fold = 1
        else:
            rep_fold = self.getValue("SELECT max(fold) FROM reproducibility WHERE Experiment = '%(track)s'")
            control_fold = self.getValue("SELECT max(fold) FROM reproducibility_vs_control")

        statement = ''' SELECT rep.Track as sample, 
                               rep.level as depth,
                               (rep.hits + 0.0)/controls.hits as ratio
                        FROM experiment_reproducibility as rep
                         INNER JOIN reproducibility_vs_control as controls
                          ON controls.Track = rep.Track AND 
                             controls.level = rep.level
                             
                        WHERE rep.Experiment = '%(track)s'  
                             AND rep.totals > 100 AND controls.fold = %(control_fold)i
                             AND rep.fold = %(rep_fold)i'''

        results = self.getAll(statement)
        results['Replicate'] = [re.match(".+\-.+\-(.+).bam", t).groups()[0] for t in results["sample"]]
        return results


class ClusterSamplesOnReproducibility(ProjectTracker):

    def __call__(self, track):

        samples = self.getDict("SELECT DISTINCT File1, totals "
                               "FROM reproducibility_distance")

        import rpy2
        ro = rpy2.robjects
        R = rpy2.robjects.r

        mat = R.matrix(ro.NA_Real, nrow=len(samples), ncol=len(samples))
        mat.rownames = ro.StrVector(samples.keys())
        mat.colnames = ro.StrVector(samples.keys())

        for sample1, sample2, hits in self.execute('''SELECT File1, File2, hits
                                                   FROM reproducibility_distance'''):
            jaccard = (hits/
                       (samples[sample1]["totals"]
                        + samples[sample2]["totals"]
                        - hits + 0.0))

            mat.rx[sample1, sample2] = jaccard

        mat_max =  R.max(mat, **{'na.rm': 'TRUE'}).ro * 1.2
        mat = mat.ro/mat_max
        distfun = R(''' function(x) as.dist(1-x) ''')

        R.x11()
        gplots = ro.packages.importr("gplots")
        gplots.heatmap_2(mat, distfun=distfun, trace="none", margins = ro.IntVector([2,10]), labCol="")
 
        return odict((("text", "#$rpl %i$#\n" % R["dev.cur"]()[0]),))


class DedupedUMIStats(ProjectTracker):

    def __call__(self,track):

        statement = ''' SELECT dus.track as sample, umi, count, (count + 0.0)/sum_count as freq
                          FROM dedup_umi_stats as dus
                       INNER JOIN (SELECT track, sum(count) as sum_count 
                                     FROM dedup_umi_stats
                                     GROUP BY track) as sum_counts
                               ON dus.track = sum_counts.track '''
        results = self.getAll(statement)
        
        results['Factor'] = [x.split("-")[0] for x in results['sample']]
        results['replicate'] = [x.split("-")[-1] for x in results['sample']]
        return results


class FragLengths(ProjectTracker):

    def getTracks(self):
        result= self.getValues(
            "SELECT DISTINCT track FROM %(table)s")
        result = [x.split(".")[0] for x in result]
        return result

    def __call__(self, track):


        statement = '''SELECT length, count
                       FROM %(table)s
                       WHERE track LIKE '%%%(track)s%%' '''


        results = self.getAll(statement)
        bins = map(int, PARAMS["experiment_length_bins"].split(","))

        bins = [bin for bin in bins if bin < max(results["Length"])]

        return_vals = OrderedDict()
        return_vals["0-%i" % bins[0]] = sum([count for start, count in
                                            zip(*results.itervalues()) if
                                           start < bins[0]])

        for i in range(len(bins)-1):
            name = "%i-%i" % (bins[i],bins[i+1])
            return_vals[name] = sum ([count for start, count in
                                      zip(*results.itervalues()) if
                                      start >= bins[i] and
                                      start < bins[i+1]])
        
        return_vals[">=%i" % bins[-1]] = sum([count for start, count in
                                              zip(*results.itervalues()) if
                                             start >= bins[-1]])

        return return_vals


class MappedFragLength(FragLengths):

    table = "mapping_frag_lengths"


class DedupedFragLengths(FragLengths):

    table = "deduped_frag_lengths"


class LengthDedupedRatios(MappedFragLength,DedupedFragLengths):

    table = "deduped_frag_lengths"

    def __call__(self, track):

        mapped = MappedFragLength()
        mapped = mapped(track)
        deduped = DedupedFragLengths()
        deduped = deduped(track)

        results = OrderedDict()

        for category in mapped:
            
            results[category] = (float(deduped[category]) /
                                 float(mapped[category]))

        return results



class TotalMergedReads(ProjectTracker):

    def getTracks(self):

        tracks = self.getValues("SELECT DISTINCT track FROM mapping.view_mapping")
        

        pattern = re.compile("merged_(.+-.+-R.+)\..+")
        tracks = [pattern.match(track).groups()[0] for track in tracks if
                  pattern.match(track)]
        
        
        return tracks


    def __call__(self,track):

        mapper = PARAMS["mappers"]
        statement = '''SELECT reads_total/2.0 as Reads
                       FROM mapping.view_mapping
                       WHERE track='merged_%(track)s.%(mapper)s' '''

        return self.getDataFrame(statement)

class SplicingIndex(iCLIPTracker):

    def getTracks(self):
        return self.getValues("SELECT track from splicing_index")

    def __call__(self, track):

        statement = '''SELECT (2*Exon_Exon+0.0)/(Exon_Intron + Intron_Exon) as SI
                       FROM splicing_index
                       WHERE track = '%(track)s' '''
        return self.getDataFrame(statement)
