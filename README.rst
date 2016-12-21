Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains:

* pipelines for re-analysis of iCLIP and scRNA-Seq data for UMI-tools publication

* ipython notebook for simulations

* ipython notebooks to generate publication figures


Single Cell Analysis
---------------------

To run the scRNA-Seq pipeline, create a new directory and copy the
files in the 'data' directory to the new directory, along with the
configuration files:

.. code:: bash

   mkdir scRNA_seq
   cp [UMI-Tool_pipelines git directory]/data/* scRNA_seq
   cp [UMI-Tool_pipelines git directory]/pipeline_scRNASeq/* scRNA_seq
   cd scRNA_seq

Then run:

.. code:: bash

   python [UMI-Tool_pipelines git directory]/pipeline_ScRNASeq.py -v10 make full

iCLIP Analysis
---------------

To run the iCLIP analysis with the default configuration, create a new
directory, download the datafiles, copy the configuration:

.. code:: bash

	  mkdir iCLIP_analysis
	  cd iCLIP_analysis
	  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR205/004/SRR2057564/SRR2057564.fastq.gz
	  .   .
	  .   .
	  .   .
	  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR205/004/SRR2057598/SRR2057598.fastq.gz
	  cp [UMI-Tools_pipelines git directory]/pipeline_iCLIP/SRSF_pipeline.ini pipeline.ini

For the complete analysis of the SRSF data download the following SRR
records from ENA: SRR2057564, SRR2057565, SRR2057566, SRR2057567,
SRR2057568, SRR2057569, SRR2057570, SRR2057571, SRR2057572,
SRR2057573, SRR20 57574, SRR2057575, SRR2057576, SRR2057577,
SRR2057578, SRR2057579, SRR2057580, SRR2057581, SRR2057582,
SRR2057583, SRR2057584, SRR2057 585, SRR2057586, SRR2057587,
SRR2057588, SRR2057589, SRR2057590, SRR2057591, SRR2057592,
SRR2057593, SRR2057594, SRR2057595, SRR2057596, SRR2057597,
SRR2057598

The pipeline requires a genome sequence (in fasta format) and a bowtie 
index, named the same way and in the same directory. The `pipeline.ini`
file must be edited to point to these files, so for example in the ini file::

	[bowtie]

	# The location of the genome and indexes for use with bowtie
	# The directory should contain a fasta genome with the .fa
	# extension and a bowtie index with the same prefix.
	# For the SRSF data set this should be an mm9 genome/index
	# fo the TDP dataset this should be a hg19 genome/index
	index_dir=~/genomes/bowtie/
	genome=mm9

Will make the pipeline look for the genome in the `~/genomes/bowtie` directory.
It will use the mm9 and expect mm9.fa and mm9.1.ebwt, mm9.2.ebwt, etc to be
present in that directory.

Now run the pipeline with:

.. code:: bash

	  python [UMI-Tool_pipelines git directory]/pipeline_iCLIP.py make full

To run the TDP analysis, copy `TDP_pipeline.ini` rather than
`SRSF_pipeline.ini`, set the bowtie index to a hg19 index and run:

.. code:: bash

	   python [UMI-Tool_pipelines git directory]/pipeline_iCLIP.py make runNotebooks1

Running locally or running on a cluster
----------------------------------------

Pipelines will run either locally or on a cluster. The default cluster
configuration is to use an SGE cluster manager with the `all.q` queue,
the `dedicated` pe environment and `mem_free`.  These defaults can be
altered by changing the following settings in the `[cluster]` section
of the pipeline.ini files:

* `queue`: The default queue to submit jobs to. Leave as NONE to have your
  cluster manager decide. (default: all.q)
* `parallel_environment`: The SGE parallel environment to request when
  submitting multi-process jobs (default: dedicated)
* `memory_resource`: The resource name to use when requesting a certian
  amount of memory for a job.  (default: mem_free)
* `pe_queue`: If this variable is set then a different queue is used
  when submitting parallel jobs. (no default)
* `options`: any other options to pass to the queue manager.

The pipelines are also compatible with the `SLURM` and `torque`
cluster managers. Set the `manager` parameter in the `[cluster]`
section of the ini. 

It is also possible to run the pipelines locally by adding
`--no-cluster` to the command.  However this will take a very long
time. Running the iCLIP pipeline on our cluster with 100 parallel jobs
(each possibly using multiple cores) takes around 50 hours. We
estimate running locally would take many weeks.

Dependencies
-------------

These pipelines require the following dependencies:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|CGAPipelines        | e6bb3be           |Pipelining infrastructure, mapping pipeline     |
|                    |                   |(http:/github.com/CGATOxford/CGATPipelines)     | 
+--------------------+-------------------+------------------------------------------------+
|CGAT                | 0.2.4             |Various                                         |
|                    |                   |(http:/github.com/CGATOxford/cgat)              |
+--------------------+-------------------+------------------------------------------------+
|Bowtie              | 1.1.2             |Mapping iCLIP reads                             |
+--------------------+-------------------+------------------------------------------------+
|BWA                 | 0.7.12-r1039      |Mapping scRNA-seq reads                         |
+--------------------+-------------------+------------------------------------------------+
|FastQC              | 0.11.2            |Quality Control of demuxed reads                |
+--------------------+-------------------+------------------------------------------------+
|bedtools            | 2.22.0            |Interval manipulation                           |
+--------------------+-------------------+------------------------------------------------+
|samtools            | 1.3.1             |Read manipulation                               |
+--------------------+-------------------+------------------------------------------------+
|UMI-tools           | 0.0.2             |UMI manipulation                                |
+--------------------+-------------------+------------------------------------------------+
|reaper              | 13-100            |Used for demuxing and clipping reads            |
+--------------------+-------------------+------------------------------------------------+
|trimmomatic         | 0.32              |Trimming reads for scRNA-seq                    |
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
|SRA toolkit         | 2.8.0             |Extracting data from SRA files                  |
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
|R                   | 3.2.1             |Figure creation (packages ggplot2, reshape,     |
|                    |                   |plyr, grid, gplots, Biobase, RColorBrewer)      |
+--------------------+-------------------+------------------------------------------------+
|jupyter             | 4.1               |Running the statistical analysis and generating |
|                    |                   |figures                                         |
+--------------------+-------------------+------------------------------------------------+
