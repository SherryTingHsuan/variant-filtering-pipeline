# Genomic Variant Filtering Pipeline

This project is an R based pipeline for filtering genomic variants in trio based sequencing data from rare disease cohorts.

The pipeline processes pre annotated input files and applies sequential filtering based on genotype, functional consequence, allele depth, population frequency, and predicted variant effect.

It also reconstructs parent child relationships from sample mapping tables and performs parent based segregation checks for homozygous and heterozygous candidate variants.

This workflow was developed during my research placement at the UCL Queen Square Institute of Neurology to improve efficiency and reproducibility in genomic data analysis.

## Key Features
- Sequential variant filtering for homozygous and heterozygous candidates
- Data cleaning and standardisation across annotated files
- Parent sample mapping using metadata tables
- Segregation checking within inferred family structure
- Structured Excel output for downstream interpretation
