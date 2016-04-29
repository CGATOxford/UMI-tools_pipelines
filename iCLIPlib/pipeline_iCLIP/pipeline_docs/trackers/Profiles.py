from iCLIPTracker import *


class TrackerImagesPlus(TrackerImages):
    
    glob=None
    pattern=None

    def __init__(self, *args, **kwargs):
        if self.glob:
            Tracker.__init__(self, *args, **kwargs)
        elif "glob" in kwargs:
            TrackerImages.__init__(self, *args, **kwargs)
        else:
            raise ValueError("Either set the self.glob or provide a :glob: parameter")

    def getPaths(self):

        filenames = glob.glob(self.glob)
        self.filenames = filenames

        if self.pattern:
            rx = re.compile(self.pattern)
            parts = [rx.search(x).groups() for x in filenames if rx.search(x)]
            self.path2file = {rx.search(x).groups(): x for x in filenames if rx.search(x)}
            result = []
            for x in range(rx.groups):
                result.append(sorted(list(set([ part[x] for part in parts ]) )))
            return result
        else:
            self.path2file = {(x,): x for x in filenames}

    def __call__(self, track, slice=None):
        
        if slice:
            key = (track, slice)
        else: 
            key = (track,)
        
        try:
            fn = self.path2file[key]
        except KeyError:
            # this particular track/slice combination doesn't exist
            return None
        
        name = " ".join(key)
        return odict((('name', name), ('filename', fn)))


class GeneProfiles(TrackerImagesPlus):
    pattern = ".+/(.+\-.+)\-(.+).tsv.geneprofilewithintrons.detail.png"


class ExonProfiles(TrackerImagesPlus):
    glob = "gene_profiles.dir/*exons.intervalprofile.detail.png"
    pattern = "./(.+\-.+)\-(.+).exons.intervalprofile.detail.png"


class IntronProfiles(TrackerImagesPlus):
    glob = "gene_profiles.dir/*introns.intervalprofile.detail.png"
    pattern = "./(.+\-.+)\-(.+).introns.intervalprofile.detail.png"


class ExonStartProfiles(TrackerImagesPlus):
    glob = "gene_profiles.dir/*exons.tssprofile.png"
    pattern = "./(.+\-.+)\-(.+).exons.tssprofile.png"


class ExonEndProfiles(TrackerImagesPlus):
    glob = "gene_profiles.dir/*introns.tssprofile.png"
    pattern = "./(.+\-.+)\-(.+).introns.tssprofile.png"


class GeneProfiles2(ProjectTracker):

    def getSlices(self):

        return self.getValues("SELECT DISTINCT rep FROM gene_profiles")

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM gene_profiles")

    def __call__(self, track, slice):

        statement = '''SELECT bin, area, region
                       FROM gene_profiles
                       WHERE rep='%(slice)s' AND factor='%(track)s' '''

        df = self.getDataFrame(statement)
        
        df["area"] = pandas.rolling_mean(df["area"], window=5)

        return df
