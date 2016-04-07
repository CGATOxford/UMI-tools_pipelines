from CGATReport.Tracker import *

import CGATPipelines.Pipeline as P

#############################################################################
# Get parameterization

P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

LOCALPARAMS = P.PARAMS

P.getParameters( ["%s/pipeline.ini" % __file__[:-len(".py")],
                  "../pipeline.ini",
                  "%s/pipeline.ini" % LOCALPARAMS["iclip_dir"],
                  "pipeline.ini" ])

PARAMS = P.PARAMS
class ProjectTracker(TrackerSQL):


    PARAMS = P.PARAMS

    def __init__(self, *args, **kwargs ):
        database_path = os.path.join(PARAMS["iclip_dir"], PARAMS["iclip_database"])
        database = database_path
        
        TrackerSQL.__init__(self, *args, backend = "sqlite:///" + database , **kwargs )
        
        # issuing the ATTACH DATABASE into the sqlalchemy ORM (self.db.execute( ... ))
        # does not work. The database is attached, but tables are not accessible in later
        # SELECT statements.

        mapping_database = os.path.join(PARAMS["iclip_dir"],"mapping.dir/csvdb")
        if not self.db:
            def _create():
                import sqlite3
                conn = sqlite3.connect(database )
                statement = '''ATTACH DATABASE '%s' as annotations;
                   ATTACH DATABASE '%s' as mapping; ''' % (self.PARAMS["annotations_database"],
                                                           mapping_database)
     
                conn.executescript(statement)
                
                return conn
              
            self.connect( creator = _create )

class iCLIPTracker(TrackerSQL):

    def __init__(self, *args, **kwargs):

        database_path = os.path.join(PARAMS["iclip_dir"], PARAMS["iclip_database"])
        mapping_database = os.path.join(PARAMS["iclip_dir"],"mapping.dir/csvdb")

        TrackerSQL.__init__(self, *args, backend = "sqlite:///" + database_path,
                            attach = ((PARAMS["annotations_database"], 'annotations'),
                                      (mapping_database, 'mapping')))
