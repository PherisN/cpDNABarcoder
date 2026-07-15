# cpDNABarcoder 1.0.1

* Fixed batch export so combined FASTA and CSV outputs retain records from all successfully processed genomes instead of overwriting earlier records.
* Added a regression test for multi-genome, multi-marker FASTA export.

# cpDNABarcoder 1.0.0

* Initial development release.
* Added chloroplast genome parsing from GenBank files.
* Added extraction support for rbcL, matK, and psbA-trnH markers.
* Added barcode QC, export, reporting, pipeline, and batch workflows.
* Added automated test coverage for core package modules.
