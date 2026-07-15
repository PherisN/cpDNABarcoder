# cpDNABarcoder

A modular R pipeline for extracting chloroplast DNA barcode loci from complete plastid genomes.

## Supported markers

- rbcL
- matK
- psbA-trnH

## Project goals

- Extract barcode loci from GenBank plastomes.
- Generate reference FASTA files.
- Produce QC reports.
- Align references with Sanger sequences.
- Support phylogenetic and DNA barcoding analyses.

## Project status

Version 1.0.1 is a stable development release that passes `R CMD check` with 0 errors, 0 warnings, and 0 notes. The package currently supports extraction of rbcL, matK, and psbA-trnH from annotated GenBank plastomes, with QC, export, reporting, and batch-processing workflows.
