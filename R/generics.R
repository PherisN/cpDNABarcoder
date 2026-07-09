###############################################################
##
## generics.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Define S3 generic accessor functions used throughout the
##   cpDNAbarcodes package.
##
## Responsibilities:
##   - Genome accessors
##   - Barcode accessors
##   - Shared object interfaces
##
###############################################################

#==============================================================
# Shared Accessors
#==============================================================

#--------------------------------------------------------------
# Retrieve object metadata
#
# Description:
#   Generic accessor that returns the metadata associated with
#   a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   A named list.
#--------------------------------------------------------------
get_metadata <- function(x){
  
  UseMethod("get_metadata")
  
}


#--------------------------------------------------------------
# Retrieve the DNA sequence
#
# Description:
#   Generic accessor that returns the DNA sequence stored in a
#   supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_sequence <- function(x){
  
  UseMethod("get_sequence")
  
}


#==============================================================
# Genome Accessors
#==============================================================

#--------------------------------------------------------------
# Retrieve the feature table
#
# Description:
#   Generic accessor that returns the parsed genomic feature
#   table associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   A data.frame.
#--------------------------------------------------------------
get_feature_table <- function(x){
  
  UseMethod("get_feature_table")
  
}

#--------------------------------------------------------------
# Retrieve the accession number
#
# Description:
#   Generic accessor that returns the GenBank accession
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_accession <- function(x){
  
  UseMethod("get_accession")
  
}

#--------------------------------------------------------------
# Retrieve the GenBank version
#
# Description:
#   Generic accessor that returns the GenBank version
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_version <- function(x){
  
  UseMethod("get_version")
  
}

#--------------------------------------------------------------
# Retrieve the genome definition
#
# Description:
#   Generic accessor that returns the GenBank definition
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_definition <- function(x){
  
  UseMethod("get_definition")
  
}

#--------------------------------------------------------------
# Retrieve the organism name
#
# Description:
#   Generic accessor that returns the scientific organism name
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_organism <- function(x){
  
  UseMethod("get_organism")
  
}

#--------------------------------------------------------------
# Retrieve the genome length
#
# Description:
#   Generic accessor that returns the total genome length
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Integer.
#--------------------------------------------------------------
get_genome_length <- function(x){
  
  UseMethod("get_genome_length")
  
}


#==============================================================
# Barcode Accessors
#==============================================================

#--------------------------------------------------------------
# Retrieve the marker name
#
# Description:
#   Generic accessor that returns the DNA barcode marker
#   represented by a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_marker <- function(x){
  
  UseMethod("get_marker")
  
}

#--------------------------------------------------------------
# Retrieve genomic coordinates
#
# Description:
#   Generic accessor that returns the genomic coordinates
#   associated with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   A named list.
#--------------------------------------------------------------
get_coordinates <- function(x){
  
  UseMethod("get_coordinates")
  
}

#--------------------------------------------------------------
# Retrieve the DNA strand
#
# Description:
#   Generic accessor that returns the DNA strand associated
#   with a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_strand <- function(x){
  
  UseMethod("get_strand")
  
}

#--------------------------------------------------------------
# Retrieve the feature type
#
# Description:
#   Generic accessor that returns the genomic feature type
#   represented by a supported object.
#
# Arguments:
#   x : Supported object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------
get_feature_type <- function(x){
  
  UseMethod("get_feature_type")
  
}