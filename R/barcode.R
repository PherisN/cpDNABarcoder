###############################################################
##
## barcode.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Construct and manage DNA barcode objects extracted from
##   chloroplast genomes.
##
## Responsibilities:
##   - Construct cpBarcode objects
##   - Validate barcode objects
##   - Provide barcode accessors
##   - Display barcode summaries
##
###############################################################

#==============================================================
# Constructors
#==============================================================

#--------------------------------------------------------------
# Create a cpBarcode object
#
# Description:
#   Constructs an S3 object representing an extracted DNA
#   barcode marker.
#
# Arguments:
#   marker       : Marker name.
#   accession    : GenBank accession.
#   version      : GenBank version.
#   organism     : Organism name.
#   feature_type : Annotation type (e.g. gene, CDS, tRNA,
#                  rRNA, intergenic).
#   sequence     : DNA sequence.
#   coordinates  : Named list containing start and end
#                  coordinates.
#   strand       : DNA strand.
#   metadata     : Optional metadata.
#
# Returns:
#   A cpBarcode object.
#--------------------------------------------------------------

new_cpBarcode <- function(marker,
                          accession,
                          version,
                          organism,
                          feature_type,
                          sequence,
                          coordinates,
                          strand,
                          metadata = list()){
  
  #---------------------------------------------
  # Validate sequence
  #---------------------------------------------
  
  if(!is_single_string(sequence)){
    
    stop(
      "sequence must be a single character string."
    )
    
  }
  
  if(nchar(sequence) == 0){
    
    stop(
      "sequence is empty."
    )
    
  }
  
  #---------------------------------------------
  # Validate coordinates
  #---------------------------------------------
  
  if(
    !is.list(coordinates) ||
    !is_valid_coordinate_pair(
      coordinates$start,
      coordinates$end
    )
  ){
    
    stop(
      "coordinates must contain 'start' and 'end'."
    )
    
  }
  
  obj <- list(
    
    marker = marker,
    
    accession = accession,
    
    version = version,
    
    organism = organism,
    
    feature_type = feature_type,
    
    sequence = sequence,
    
    coordinates = coordinates,
    
    strand = strand,
    
    length = nchar(sequence),
    
    metadata = metadata
    
  )
  
  class(obj) <- "cpBarcode"
  
  return(
    validate_barcode(obj)
  )
  
}


#--------------------------------------------------------------
# Validate a cpBarcode input object
#
# Description:
#   Ensures that an object inherits from the cpBarcode class.
#
# Arguments:
#   barcode : Object to validate.
#
# Returns:
#   Invisibly returns the validated object.
#--------------------------------------------------------------

assert_cpBarcode <- function(barcode){
  
  if(!inherits(barcode, "cpBarcode")){
    
    stop(
      "Input must be a cpBarcode object."
    )
    
  }
  
  invisible(barcode)
  
}


#==============================================================
# Accessors
#==============================================================

#--------------------------------------------------------------
# Retrieve the DNA sequence from a cpBarcode object
#
# Description:
#   Returns the nucleotide sequence stored in a cpBarcode object.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string containing the DNA sequence.
#--------------------------------------------------------------

get_sequence.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$sequence)
  
}


#--------------------------------------------------------------
# Retrieve genomic coordinates from a cpBarcode object
#
# Description:
#   Returns the genomic coordinates associated with a barcode.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   A named list containing the start and end coordinates.
#--------------------------------------------------------------

get_coordinates.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$coordinates)
  
}


#--------------------------------------------------------------
# Retrieve metadata from a cpBarcode object
#
# Description:
#   Returns the metadata associated with a cpBarcode object.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   A named list.
#--------------------------------------------------------------

get_metadata.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$metadata)
  
}


#--------------------------------------------------------------
# Retrieve the marker name from a cpBarcode object
#
# Description:
#   Returns the DNA barcode marker represented by the object.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_marker.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$marker)
  
}


#--------------------------------------------------------------
# Retrieve the accession number from a cpBarcode object
#
# Description:
#   Returns the GenBank accession associated with a barcode.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_accession.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$accession)
  
}


#--------------------------------------------------------------
# Retrieve the organism name from a cpBarcode object
#
# Description:
#   Returns the scientific name of the organism associated with
#   a barcode.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_organism.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$organism)
  
}


#--------------------------------------------------------------
# Retrieve the GenBank version from a cpBarcode object
#
# Description:
#   Returns the GenBank version associated with a barcode.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_version.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$version)
  
}


#--------------------------------------------------------------
# Retrieve the feature type from a cpBarcode object
#
# Description:
#   Returns the genomic feature type represented by the
#   barcode (for example gene, CDS, tRNA, rRNA, or
#   intergenic).
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_feature_type.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$feature_type)
  
}


#--------------------------------------------------------------
# Retrieve the strand from a cpBarcode object
#
# Description:
#   Returns the DNA strand on which the barcode is located.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   Character string ("+" or "-").
#--------------------------------------------------------------

get_strand.cpBarcode <- function(x){
  
  assert_cpBarcode(x)
  
  return(x$strand)
  
}

#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate a cpBarcode object
#
# Description:
#   Ensures that a cpBarcode object contains all required fields
#   and a non-empty DNA sequence.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   The validated cpBarcode object.
#--------------------------------------------------------------

validate_barcode <- function(barcode){
  
  #---------------------------------------------
  # Validate object class
  #---------------------------------------------
  
  if(!inherits(barcode, "cpBarcode")){
    
    stop(
      "Input must be a cpBarcode object."
    )
    
  }
  
  #---------------------------------------------
  # Validate coordinates
  #---------------------------------------------
  
  if(
    !is.list(barcode$coordinates) ||
    !is_valid_coordinate_pair(
      barcode$coordinates$start,
      barcode$coordinates$end
    )
  ){
    
    stop(
      "Barcode coordinates are invalid."
    )
    
  }
  
  if(
    barcode$coordinates$end <
    barcode$coordinates$start
  ){
    
    stop(
      "Barcode coordinates are inconsistent."
    )
    
  }
  
  #---------------------------------------------
  # Validate required fields
  #---------------------------------------------
  
  required <- c(
    
    "marker",
    "accession",
    "version",
    "organism",
    "feature_type",
    "sequence",
    "coordinates",
    "strand",
    "length",
    "metadata"
    
  )
  
  missing <- required[
    !required %in% names(barcode)
  ]
  
  if(length(missing) > 0){
    
    stop(
      "Missing barcode fields: ",
      paste(missing, collapse = ", ")
    )
    
  }
  
  #---------------------------------------------
  # Validate sequence
  #---------------------------------------------
  
  if(
    !is_single_string(barcode$sequence) ||
    !is_valid_dna(barcode$sequence)
  ){
    
    stop(
      "Barcode sequence is empty."
    )
    
  }
  
  if(barcode$length != nchar(barcode$sequence)){
    
    stop(
      "Stored barcode length does not match the sequence length."
    )
    
  }
  
  #---------------------------------------------
  # Validate strand
  #---------------------------------------------
  
  if(!(barcode$strand %in% c(
    STRAND_FORWARD,
    STRAND_REVERSE
  ))) {
    
    stop(
      "Barcode strand is invalid."
    )
    
  }
  
  #---------------------------------------------
  # Validate feature type
  #---------------------------------------------
  
  valid_types <- c(
    
    FEATURE_GENE,
    
    FEATURE_CDS,
    
    FEATURE_TRNA,
    
    FEATURE_RRNA,
    
    FEATURE_INTERGENIC
    
  )
  
  if(!(barcode$feature_type %in% valid_types)){
    
    warning(
      "Unknown feature type."
    )
    
  }
  
  #---------------------------------------------
  # Validate marker
  #---------------------------------------------
  
  valid_markers <- vapply(
    
    MARKER_DICTIONARY,
    
    function(x) x$marker,
    
    character(1)
    
  )
  
  if(!(barcode$marker %in% valid_markers)){
    
    warning(
      "Unknown barcode marker."
    )
    
  }
  
  #---------------------------------------------
  # Validate metadata
  #---------------------------------------------
  
  if(!is.list(barcode$metadata)){
    
    stop(
      "Metadata must be a list."
    )
    
  }
  
  return(barcode)
  
}


#==============================================================
# Methods
#==============================================================

#--------------------------------------------------------------
# Print a cpBarcode object
#
# Description:
#   Displays a concise summary of a cpBarcode object.
#
# Arguments:
#   x : A cpBarcode object.
#
# Returns:
#   Invisibly returns the object.
#--------------------------------------------------------------

print.cpBarcode <- function(x, ...){
  
  x <- validate_barcode(x)
  
  cat("\n")
  
  cat("cpBarcode Object\n")
  
  cat("---------------------------\n")
  
  cat("Marker       :", x$marker, "\n")
  
  cat("Feature type :", x$feature_type, "\n")
  
  cat("Organism     :", x$organism, "\n")
  
  cat("Accession    :", x$accession, "\n")
  
  cat("Version      :", x$version, "\n")
  
  cat(
    "Coordinates  :",
    paste(
      x$coordinates$start,
      x$coordinates$end,
      sep = "-"
    ),
    "\n"
  )
  
  cat("Length       :", x$length, "bp\n")
  
  cat("Strand       :", x$strand, "\n")
  
  invisible(x)
  
}