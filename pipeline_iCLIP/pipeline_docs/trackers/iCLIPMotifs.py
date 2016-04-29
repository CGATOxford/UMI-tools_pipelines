from iCLIPTracker import *
import xml.etree.ElementTree
import re

def re2concensus(motif):
    motif = re.sub("\n", "", motif)
    pattern = "\[(.)[^\]]+\]"
    return re.sub(pattern, r"\1", motif)

class MemeResults(ProjectTracker):

    trees = {}

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM meme_summary")

    def getSlices(self):

        tracks = self.getTracks()

        slices = []

        for track in tracks:
            
            resultsdir = os.path.abspath(
                os.path.join(PARAMS["iclip_exportdir"], "meme", "%s.meme" % track))
            if not os.path.exists(resultsdir):
                continue

            tree = xml.etree.ElementTree.ElementTree()
            tree.parse(os.path.join(resultsdir, "meme.xml"))
            self.trees[track] = tree

            motifs = tree.find("motifs")
            for motif in motifs.getiterator("motif"):
                regex = motif.find("regular_expression").text
                slices.append(re2concensus(regex))

        return list(set(slices))

    def __call__(self, track, slice=None):

        tree = self.trees.get(track, None)
        if tree is None:
            return []
        resultsdir = os.path.abspath(
                os.path.join(PARAMS["iclip_exportdir"], "meme", "%s.meme" % track))

        sequences = tree.find("model").find("num_sequences").text
        motifs = tree.find("motifs")
        nmotif = 0
        result = odict()
        for motif in motifs.getiterator("motif"):
            
            nmotif += 1
            regex = motif.find("regular_expression").text
            if not re2concensus(regex) == slice:
                continue
            motif_img = "%s/logo%i.png" % (resultsdir, nmotif)
            motif_rc_img = "%s/logo_rc%i.png" % (resultsdir, nmotif)
            img, rc_img = "na", "na"
            if os.path.exists(motif_img):
                img = '''.. image:: %s
   :scale: 25%% ''' % motif_img
            if os.path.exists(motif_rc_img):
                rc_img = '''.. image:: %s
   :scale: 25%% '''  % motif_rc_img

            result[str(nmotif)] = odict((
                ("width", motif.get("width")),
                ("evalue", motif.get("e_value")),
                ("information content", motif.get("ic")),
                ("sites", "%s/%s" % (motif.get("sites"), sequences)),
                ("link", "`meme_%s_%i <%s/meme.html#summary%i>`_" %
                 (track, nmotif, resultsdir, nmotif)),
                ("img", img),
               
            ))

        if len(result) > 0:
            return result

class DremeResults(ProjectTracker):

    trees = {}

    def getTracks(self):
        tracks = glob.glob(os.path.abspath(os.path.join(PARAMS["iclip_dir"],
                                 self.glob_files)))
        print "glob returned %i tracks" % len(tracks)
        tracks = [re.search(self.glob_pattern, track).groups()[0] for track in tracks]

        return tracks

    def getSlices(self):

        tracks = self.getTracks()
        slices = []

        for track in tracks:
            pattern = re.sub("\(\.\+\)", "%s", self.glob_pattern)

            resultsdir = os.path.abspath(
                os.path.join(PARAMS["iclip_dir"], pattern % track))
            if not os.path.exists(resultsdir):
                print "file %s does not exist" % resultsdir
                continue

            tree = xml.etree.ElementTree.ElementTree()
            tree.parse(os.path.join(resultsdir, "dreme.xml"))
            print "adding %s to trees" % track
            self.trees[track] = tree

            motifs = tree.find("motifs")
            for motif in motifs.getiterator("motif"):
                slices.append(motif.get("seq"))

        return list(set(slices))

    def __call__(self, track, slice=None):

        pattern = re.sub("\(\.\+\)", "%s", self.glob_pattern)
        resultsdir = os.path.abspath(os.path.join(PARAMS["iclip_dir"],
            os.path.join(pattern % track)))

        if not os.path.exists(resultsdir):
            print "resultsdir %s does not exist"
            return []

        tree = self.trees[track]

        if tree is None:
            return []

        model = tree.find("model")
        num_positives = int(model.find("positives").get("count"))
        num_negatives = int(model.find("negatives").get("count"))

        motifs = tree.find("motifs")
        nmotif = 0
        result = pandas.DataFrame(columns = ["sequence","evalue", "positives",
                                             "negatives","enrichment", "link",
                                             "img"])

        for motif in motifs.getiterator("motif"):

            seq = motif.get("seq")
            if not seq == slice:
                continue

            nmotif += 1
            id = motif.get("id")
            seq = motif.get("seq")

            motif_img = "%(resultsdir)s/%(id)snc_%(seq)s.png" % locals()
            img, rc_img = "na", "na"
            if os.path.exists(motif_img):
                img = '''.. image:: %s
   :scale: 25%% ''' % motif_img

            p = float(motif.get("p"))
            n = float(motif.get("n"))
            
            try:
                enrichment = (p/num_positives)/(n/num_negatives)
                if enrichment < self.enrichment_threshold:
                    continue
                enrichment = "{:.0%}".format(
                    enrichment)
            except ZeroDivisionError:
                enrichment = "inf"

            result = result.append(dict((
                ("sequence", seq),
                ("evalue", motif.get("evalue")),
                ("positives", "{:d}/{} ({:.0%})".format(int(p), num_positives,
                                                        p/num_positives)),
                ("negatives", "{:d}/{} ({:.0%})".format(int(n), num_negatives,
                                                        n/num_negatives)),
                ("enrichment", enrichment),
                ("link", "`dreme_%s <%s/dreme.html>`_" %
                 (track, resultsdir)),
                ("img", img),
            )), ignore_index=True)

        if result.shape[0] > 0:
            return result


class SimpleDremeResults(DremeResults):
    enrichment_threshold = 1.5
    glob_files = "export/dreme/*"
    glob_pattern = "export/dreme/(.+)"
