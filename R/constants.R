###############################################################
##
## constants.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Define package-wide constants used throughout
##   cpDNAbarcodes.
##
## Responsibilities:
##   - GenBank record identifiers
##   - Feature type constants
##   - Marker class constants
##   - Strand constants
##   - DNA alphabet constants
##   - Export defaults
##
###############################################################


#==============================================================
# Package Information
#==============================================================

PACKAGE_NAME <- "cpDNAbarcodes"

PACKAGE_VERSION <- "1.0.0"


#==============================================================
# GenBank Section Headers
#==============================================================

LOCUS_START <- "^LOCUS"

DEFINITION_START <- "^DEFINITION"

ACCESSION_START <- "^ACCESSION"

VERSION_START <- "^VERSION"

ORGANISM_START <- "^  ORGANISM"

GENBANK_FEATURE_START <- "^FEATURES"

ORIGIN_START <- "^ORIGIN"

GENBANK_END <- "^//"


#==============================================================
# Feature Types
#==============================================================

FEATURE_GENE <- "gene"

FEATURE_CDS <- "CDS"

FEATURE_TRNA <- "tRNA"

FEATURE_RRNA <- "rRNA"

FEATURE_INTERGENIC <- "intergenic"


#==============================================================
# Marker Classes
#==============================================================

MARKER_CLASS_GENE <- "GENE"

MARKER_CLASS_INTERGENIC <- "INTERGENIC"

MARKER_CLASS_NUCLEAR <- "NUCLEAR"


#==============================================================
# Strand Constants
#==============================================================

STRAND_FORWARD <- "+"

STRAND_REVERSE <- "-"


#==============================================================
# Feature Location Types
#==============================================================

LOCATION_SIMPLE <- "simple"

LOCATION_JOIN <- "join"

LOCATION_ORDER <- "order"


#==============================================================
# DNA Alphabet
#==============================================================

DNA_BASES <- c(
  "A",
  "C",
  "G",
  "T"
)

IUPAC_BASES <- c(
  DNA_BASES,
  "R", "Y", "S", "W",
  "K", "M", "B", "D",
  "H", "V", "N"
)


#==============================================================
# Export Defaults
#==============================================================

DEFAULT_FASTA_WIDTH <- 80

DEFAULT_FASTA_EXTENSION <- ".fasta"
