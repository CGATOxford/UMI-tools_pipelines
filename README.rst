Tools for dealing with Unique Molecular Identifiers
====================================================
This repository contains:

* pipelines for re-analysis of iCLIP and scRNA-Seq data for UMI-tools publication

* ipython notebook for simulations

* ipython notebooks to generate publication figures


To run the scRNA-Seq pipeline, create a new directory and copy the files in the 'data' directory to the new directory, along with the configuration files:

.. code:: bash

   mkdir scRNA_seq
   cp [UMI-Tool_pipelines git directory]/data/* scRNA_seq
   cp [UMI-Tool_pipelines git directory]/pipeline_scRNASeq/* scRNA_seq
   cd scRNA_seq

Then run:

.. code:: bash

   python [UMI-Tool_pipelines git directory]/PipelineScRNASeq.py -v10 make full

