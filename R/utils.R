###############################################################
##
## utils.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Provide general-purpose utility functions used throughout
##   the cpDNAbarcodes package.
##
## Responsibilities:
##   - String utilities
##   - Coordinate utilities
##   - DNA sequence utilities
##   - File utilities
##
###############################################################

#==============================================================
# String Utilities
#==============================================================

#--------------------------------------------------------------
# Check whether an object is a single character string
#
# Description:
#   Returns TRUE if the supplied object is a non-missing
#   character vector of length one.
#
# Arguments:
#   x : Object to test.
#
# Returns:
#   Logical.
#--------------------------------------------------------------

is_single_string <- function(x){
  
  is.character(x) &&
    length(x) == 1 &&
    !is.na(x)
  
}

#--------------------------------------------------------------
# Check whether an object is a non-empty character vector
#
# Description:
#   Returns TRUE if the supplied object is a character vector
#   containing one or more elements.
#
# Arguments:
#   x : Object to test.
#
# Returns:
#   Logical.
#--------------------------------------------------------------

is_character_vector <- function(x){
  
  is.character(x) &&
    length(x) > 0
  
}

#--------------------------------------------------------------
# Convert a string to uppercase after trimming whitespace
#
# Description:
#   Removes leading and trailing whitespace and converts the
#   resulting string to uppercase.
#
# Arguments:
#   x : A non-missing character string of length one.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

trim_and_upper <- function(x){
  
  if(!is_single_string(x)){
    
    stop(
      "Input must be a non-missing character string of length one."
    )
    
  }
  
  toupper(
    trimws(x)
  )
  
}

#==============================================================
# Coordinate Utilities
#==============================================================

#--------------------------------------------------------------
# Check whether a coordinate is valid
#
# Description:
#   Tests whether a genomic coordinate is a positive integer.
#
# Arguments:
#   x : Coordinate.
#
# Returns:
#   Logical.
#--------------------------------------------------------------

is_valid_coordinate <- function(x){
  
  is.numeric(x) &&
    length(x) == 1 &&
    !is.na(x) &&
    x >= 1 &&
    x == as.integer(x)
  
}


#--------------------------------------------------------------
# Validate a coordinate pair
#
# Description:
#   Checks that a start and end coordinate are valid.
#
# Arguments:
#   start : Start coordinate.
#   end   : End coordinate.
#
# Returns:
#   TRUE if valid.
#--------------------------------------------------------------

is_valid_coordinate_pair <- function(start,
                                     end){
  
  is_valid_coordinate(start) &&
    is_valid_coordinate(end)
  
}


#==============================================================
# DNA Utilities
#==============================================================

#--------------------------------------------------------------
# Compute the DNA complement
#
# Description:
#   Returns the nucleotide complement of a DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

complement_sequence <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  chartr(
    
    "ACGTRYKMSWBDHVN",
    
    "TGCAYRMKSWVHDBN",
    
    toupper(sequence)
    
  )
  
}


#--------------------------------------------------------------
# Compute the reverse complement
#
# Description:
#   Returns the reverse complement of a DNA sequence.
#   Reverse complementation is only performed when the sequence
#   originates from the reverse strand.
#
# Arguments:
#   sequence : DNA sequence.
#   strand   : "+" or "-".
#
# Returns:
#   Character string.
#--------------------------------------------------------------

reverse_complement <- function(sequence){
  
  complement <- complement_sequence(sequence)
  
  bases <- strsplit(
    complement,
    ""
  )[[1]]
  
  paste(
    rev(bases),
    collapse = ""
  )
  
}


#--------------------------------------------------------------
# Return DNA sequence length
#
# Description:
#   Returns the length of a validated DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Integer.
#--------------------------------------------------------------

dna_length <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  nchar(sequence)
  
}


#==============================================================
# File Utilities
#==============================================================

#--------------------------------------------------------------
# Ensure that an output directory exists
#
# Description:
#   Creates a directory (including parent directories) if it
#   does not already exist.
#
# Arguments:
#   path : Directory path.
#
# Returns:
#   Invisibly returns the directory path.
#--------------------------------------------------------------

ensure_directory <- function(path){
  
  if(!dir.exists(path)){
    
    success <- dir.create(
      path,
      recursive = TRUE
    )
    
    if(!success){
      
      stop(
        "Unable to create directory: ",
        path
      )
      
    }
    
  }
  
  invisible(path)
  
}