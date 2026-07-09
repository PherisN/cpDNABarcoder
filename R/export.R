###############################################################
##
## export.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Export DNA barcode objects to standard bioinformatics
##   formats.
##
## Responsibilities:
##   - Format FASTA output
##   - Export FASTA files
##   - Convert barcode objects to tabular data
##   - Export barcode objects
##
###############################################################

#==============================================================
# FASTA Utilities
#==============================================================

#--------------------------------------------------------------
# Wrap a DNA sequence to a fixed line width
#
# Description:
#   Splits a DNA sequence into fixed-width lines suitable for
#   FASTA output.
#
# Arguments:
#   sequence : DNA sequence.
#   width    : Maximum number of bases per line.
#
# Returns:
#   Character vector.
#--------------------------------------------------------------

wrap_fasta_sequence <- function(sequence,
                                width = 80){
  
  if(!is_valid_dna(sequence)){
    stop(
      "Sequence contains invalid DNA symbols."
    )
  }
  
  if(dna_length(sequence) == 0){
    
    stop(
      "Sequence is empty."
    )
    
  }
  
  if(width <= 0){
    
    stop(
      "Width must be greater than zero."
    )
    
  }
  
  n <- nchar(sequence)
  
  starts <- seq(
    1,
    n,
    by = width
  )
  
  ends <- pmin(
    starts + width - 1,
    n
  )
  
  wrapped <- mapply(
    
    substr,
    
    x = sequence,
    
    start = starts,
    
    stop = ends,
    
    USE.NAMES = FALSE
    
  )
  
  wrapped
  
}

#--------------------------------------------------------------
# Build a FASTA header
#
# Description:
#   Constructs a standardized FASTA header for a cpBarcode
#   object.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string beginning with ">".
#--------------------------------------------------------------

build_fasta_header <- function(barcode){
  
  barcode <- validate_barcode(barcode)
  
  accession <- get_accession(barcode)
  
  marker <- get_marker(barcode)
  
  organism <- get_organism(barcode)
  
  coords <- get_coordinates(barcode)
  
  paste(
    
    ">",
    
    accession,
    
    "|",
    
    marker,
    
    "|",
    
    organism,
    
    "|",
    
    coords$start,
    
    "-",
    
    coords$end,
    
    sep = ""
    
  )
  
}

#--------------------------------------------------------------
# Convert a barcode to FASTA
#
# Description:
#   Converts a cpBarcode object into FASTA-formatted text.
#
# Arguments:
#   barcode : A cpBarcode object.
#   width   : FASTA line width.
#
# Returns:
#   Character vector.
#--------------------------------------------------------------

as_fasta <- function(barcode,
                     width = DEFAULT_FASTA_WIDTH){
  
  barcode <- validate_barcode(barcode)
  
  header <- build_fasta_header(barcode)
  
  sequence <- wrap_fasta_sequence(
    
    get_sequence(barcode),
    
    width
    
  )
  
  c(
    header,
    sequence
  )
  
}

#==============================================================
# FASTA Export
#==============================================================

#--------------------------------------------------------------
# Write a FASTA file
#
# Description:
#   Writes a cpBarcode object to a FASTA file.
#
# Arguments:
#   barcode : A cpBarcode object.
#   file    : Output filename.
#   width   : FASTA line width.
#
# Returns:
#   Invisibly returns the output filename.
#--------------------------------------------------------------

write_fasta <- function(barcode,
                        file,
                        width = 80){
  
  barcode <- validate_barcode(barcode)
  
  directory <- dirname(file)
  
  ensure_directory(
    dirname(file)
  )
  
  writeLines(
    
    as_fasta(
      barcode,
      width
    ),
    
    con = file
    
  )
  
  invisible(file)
  
}

#==============================================================
# Tabular Export
#==============================================================

#--------------------------------------------------------------
# Convert a barcode to a data frame
#
# Description:
#   Converts a cpBarcode object into a one-row data frame.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   data.frame.
#--------------------------------------------------------------

as_dataframe <- function(barcode){
  
  barcode <- validate_barcode(barcode)
  
  coords <- get_coordinates(barcode)
  
  data.frame(
    
    marker = get_marker(barcode),
    
    accession = get_accession(barcode),
    
    version = get_version(barcode),
    
    organism = get_organism(barcode),
    
    feature_type = get_feature_type(barcode),
    
    start = coords$start,
    
    end = coords$end,
    
    strand = get_strand(barcode),
    
    length = dna_length(get_sequence(barcode)),
    
    sequence = get_sequence(barcode),
    
    stringsAsFactors = FALSE
    
  )
  
}

#--------------------------------------------------------------
# Write a CSV file
#
# Description:
#   Writes a cpBarcode object to a CSV file.
#
# Arguments:
#   barcode : A cpBarcode object.
#   file    : Output filename.
#
# Returns:
#   Invisibly returns the output filename.
#--------------------------------------------------------------

write_csv <- function(barcode,
                      file){
  
  barcode <- validate_barcode(barcode)
  
  directory <- dirname(file)
  
  ensure_directory(
    dirname(file)
  )
  
  write.csv(
    
    as_dataframe(barcode),
    
    file,
    
    row.names = FALSE
    
  )
  
  invisible(file)
  
}

#==============================================================
# High-level API
#==============================================================

#--------------------------------------------------------------
# Write a barcode object
#
# Description:
#   Writes a cpBarcode object using the requested output
#   format.
#
# Arguments:
#   barcode : A cpBarcode object.
#   file    : Output filename.
#   format  : Output format.
#   ...     : Additional arguments passed to the exporter.
#
# Returns:
#   Invisibly returns the output filename.
#--------------------------------------------------------------

write_barcode <- function(barcode,
                          file,
                          format = "fasta",
                          ...){
  
  barcode <- validate_barcode(barcode)
  
  format <- tolower(
    trimws(format)
  )
  
  format <- sub(
    "^\\.",
    "",
    format
  )
  
  result <- switch(
    
    format,
    
    fasta = write_fasta(
      barcode,
      file,
      ...
    ),
    
    csv = write_csv(
      barcode,
      file
    ),
    
    stop(
      "Unsupported export format: ",
      format
    )
    
  )
  
  invisible(result)
  
}