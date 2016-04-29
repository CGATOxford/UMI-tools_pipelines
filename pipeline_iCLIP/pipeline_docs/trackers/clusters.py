from iCLIPTracker import *
from statsmodels.stats.multitest import multipletests
from Sample_QC import ContextStats

class ClusterStats(TrackerDataframes):

    def __call__(self, track):
        print "Reading %s" % track
        data = pandas.read_csv(self.openFile(track),
                               header=0,
                               names=["contig", "start", "p"],
                               sep="\t")
        print "Done"
        data["qvalues"] = multipletests(data["p"], method="fdr_bh")[1]

        output = dict()

        output["Bases"] = data.shape[0]
        output["Significant"] = (data["qvalues"] < 0.01).sum()
        output["Fraction_Significant"] = \
            float(output["Significant"])/output["Bases"]

        return output


class ClusterCounts(TrackerSQL):

    def __call__(self, track):

        statement = '''SELECT sample, replicate, count FROM cluster_counts'''

        return self.getDataFrame(statement)


class ClusterContexts(ContextStats):
    method="cluster"

