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

Version 1.0.1 is complete and passes `R CMD check` with 0 errors, 0 warnings, and 0 notes. This patch release fixes combined batch export so multi-genome marker FASTA and CSV outputs retain all successfully processed records.
