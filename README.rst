.. |date| date::

******************
STAR-RSEM pipeline
******************

:authors: Komal Rathi
:contact: rathik@email.chop.edu
:organization: DBHi, CHOP
:status: Completed
:date: |date|

.. meta::
   :keywords: star, rsem, 2016
   :description: DBHi STAR-RSEM pipeline.

Introduction
============

STAR-RSEM pipeline based on: 
https://github.com/BD2KGenomics/toil-scripts/blob/releases/2.0.x/src/toil_scripts/rnaseq_cgl/README.md

Annotations, Reference genome and Pipeline used:
https://github.com/BD2KGenomics/toil-rnaseq

This pipeline was used by UCSC to generate STAR-RSEM normalized counts for big datasets like TCGA, GTEx, PNOC and TARGET. (bioRxiv_)

There are five steps in this pipeline:

1. STAR genome generation
2. RSEM genome generation
3. SRA to FASTQ conversion
4. STAR alignment
5. RSEM normalization. This step creates gene and transcript level quantifications.

Installation
============

Create conda environment:

.. code-block:: bash

	conda create --name rnaseq-env
	source activate rnaseq-env
	conda install -c biobuilds sra-tools=2.5.6
	conda install -c bioconda rsem=1.2.28
	conda install -c bioconda star=2.5.2b

How to run the pipeline
=======================

.. code-block:: bash

	1. First change the **reference freeze** in the config file depending on your dataset.
	2. Run on a SGE cluster (or modify according to your system)
	
	snakemake -p -j 10 -s Snakefile --configfile config.yaml --cluster-config cluster.yaml -c "qsub -cwd -e error.txt -o output.txt -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}" &

.. links:
.. _bioRxiv: https://www.biorxiv.org/content/early/2016/07/07/062497