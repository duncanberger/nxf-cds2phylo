# nxf-cds2phylo [Pre-release]
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.4-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
## Table of contents 
* [Introduction](#Introduction)
* [Pipeline summary](#pipeline_summary)
* [Installation](#install)
* [Running cds2phylo](#run)


## Introduction <a name="Introduction"></a>
**cds2phylo** is a Nextflow pipeline for to create phylogenies based on alignments of large number of (core) genes. 

**[Note: cds2phylo is currently being developed and has only been run internally]** <br />

**[Note: Addition of executors to the pipeline is in progress - most likely will default to Slurm]**

## Pipeline summary <a name="pipeline_summary"></a>

The pipeline takes a FASTA format file of unaligned nucleotide sequences (CDSs) for each gene/allele family and each sample. Then for every gene family, it aligns the CDS sequences across all sample where the sequence is provided. The alignments are trimmed and the remaining sequence are eiether concatenated into a supermatrix or retained as individual alignments (e.g. for partitioned analysis). Finally, phylogenetic inference is performed using RAXML, IQTREE or FastTree, depending on user specifications. 

## Installation <a name="install"></a>
#### Clone repository
```
git clone https://github.com/duncanberger/nxf-cds2phylo.git
```
#### Create a virtual environment
```
mamba env create -f env.yaml
```
## Running cds2phylo <a name="run"></a>

#### Usage
```
Usage:
        nextflow run cds2phylo.nf  --fasta <FASTA> --input <GENE LIST>

Mandatory arguments:
        --input                         Path to file containing list of geneIDs to include
        --fasta                         Path to file containing fasta sequences to process

Optional arguments:
        --prefix                        Output file prefix ["out"]
        --outdir                        Output directory ["results"]
        --phylo_model                   Use Fasttree or IQTREE ["fasttree"]
        --iqtree_parameters             Other parameters to pass to IQTREE [""]
        --partition                     Do partitioned analysis instead of supermatrix [false]
        --iqtree_model_supermatrix      IQTREE supermatrix analysis model ["MFP+ASC"]
        --iqtree_model_partition        IQTREE partioned analysis model ["MFP"]
        --skip_tree                     Skips creating the phylogeny [false]
```
#### Single phylogeny
```
nextflow run cds2phylo.nf --phylo_method iqtree --fasta <FASTA> --input <GENE LIST>
```
#### One phylogeny per-sample
```
nextflow run drak.nf --prefix <SAMPLE_ID> --phylo_method iqtree --fasta <FASTA> --input <GENE LIST>
```
