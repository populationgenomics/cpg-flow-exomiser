# cpg-flow exomiser

An extraction and reimplementation of the Exomiser workflow implemented at The Centre for Population Genomics.

Original source in Production-Pipelines [Here](https://github.com/populationgenomics/production-pipelines/blob/b79e5b5c6550db138244094f89e6e244972540da/cpg_workflows/stages/exomiser.py)

This repository contains a workflow built using the [CPG-Flow framework](https://github.com/populationgenomics/cpg-flow), which is a Hail Batch based workflow executor.

## Purpose

This workflow generates Exomiser results for families.

It does this by identifying families which have not previously received an analysis, creating a pseudo-joint-call from all members of the family, pulling all the family members' metadata into a phenopacket, and running Exomiser.

Multiple single-family analyses are stacked into a single VM and executed in parallel, to make best use of the localised annotation resources (large files). 

Once analysis is complete, the resulting Gene- and Variant-level result TSVs are aggregated. 

* The Gene-level results are placed into a Seqr-ready format to be used in enhanced gene searches
* The Variant-level results are placed into a Talos-ready format to be used in variant prioritisation

## Structure

```commandline
src
├── cpg_exomiser
│   ├── __init__.py
│   ├── config_template.toml
│   ├── jobs
│   │   ├── CreateFamilyVcfs.py
│   │   ├── MakePedExtracts.py
│   │   ├── MakePhenopackets.py
│   │   └── RunExomiser.py
│   ├── run_workflow.py
│   ├── scripts
│   │   ├── combine_exomiser_gene_tsvs.py
│   │   └── combine_exomiser_variant_tsvs.py
│   ├── stages.py
```
