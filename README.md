# BCB570 project

## Instructions on running this script

Every step of this workflow is coordinated with snakemake. You can install
snakemake with pip as follows.

```
pip3 install snakemake
```

## Introduction


Do the orphan genes cluster in their own little networks? If so, are these
functional clusters or artifacts of common origin? For example, were they
all animated by insertions by the same type of telomere?

1. Arabidopsis thaliana RNA-seq: take maybe 3 studies. Run through my
   RNA-seq pipeline using Kallisto. Process with Sleuth. Include the
   bootstrapping procedure.

1. Use snakefile to organize run

1. Compare to interaction network data, using IntAct

1. Compare to promoter motif set similarity. How to do this? Retrieve
   putative set of discrete regulatory elements (one exists for
   Arabidopsis), Draw edges between genes if there share a common element.
   Weights are relative to the number of shared elements. This is not the
   only way to do it, but it may work.

1. Compare genomic distance network.

Questions:

1. Do the orphans cluster together in expression?

1. If so, do they cluster more than I would expect?

1. Take the largest or most interesting cluster

    1. Do they share more interactions than I would expect? If so, this is
       evidence that the expression clustering has some basis in protein
       function, rather than chance.
    
    1. Do they share a common, simple promoter motif?

    1. Do they share a common upstream promoter?

Naive hypothesis:

If the genes are co-expressed they must be part of the same functional
pathway. Expression is driven by selection.

Other possibilities:

They share a common origin. Expression is neutral. The co-expression
predates any shared functionality.

Story 1: A transposon proliferation provided the promoter machinery that
animated the non-genic precursor of the orphan.

Story 2: If the promoter is completely de novo, the most likely regulatory
motifs to have arisen are the shortest, most simple ones.

Perhaps, because they are coexpressed originally, and since they therefore
have the chance to interact, selection will work to enhance the
similarities, thus drawing the cluster together.


Orphan problems:

    Generally GO is not an option. There are only a few orphans that even have
    GO terms.

    Generally, evolutionary methods are not available.

    There is a few added layer of uncertainty: the orphans may not actually
    express a functional (whatever that means) protein product. They may be
    noise.


# Selection of runs
 
From the selected studies, I removed any runs that

 1. Where not contentional transcriptomic (for example, bisulfite sequencing)
 1. Where not from Col-0 *A. thaliana*

# Data

bindingsite-data.tbl - experimentally confirmed [TF-binding
sites](http://arabidopsis.med.ohio-state.edu/AtcisDB/bindingsites.html) from
Agris.

orphans - an in-house list of orphan TAIR10 model ids.

at.fa - CDS for arabidopsis models (retrieved from the Kallisto website, and
they got it from Ensembl)
