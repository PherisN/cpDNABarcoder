###############################################################
##
## io.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Read GenBank flat files and extract major sections used
##   throughout the cpDNAbarcodes pipeline.
##
## Responsibilities:
##   - Read GenBank records
##   - Extract accession metadata
##   - Extract organism metadata
##   - Extract FEATURES section
##   - Extract genome sequence
##
###############################################################

#--------------------------------------------------------------
# Read a GenBank file
#
# Description:
#   Reads a GenBank flat file into a character vector.
#
# Arguments:
#   file : Path to a GenBank file.
#
# Returns:
#   Character vector containing the GenBank record.
#--------------------------------------------------------------

read_genbank <- function(file){
  
  if(!is_single_string(file)){
    
    stop(
      "'file' must be a non-missing character string."
    )
    
  }
  
  if(!file.exists(file)){
    
    stop(
      "GenBank file not found: ",
      file
    )
    
  }
  
  if(file.info(file)$isdir){
    
    stop(
      "'file' must refer to a GenBank file, not a directory."
    )
    
  }
  
  readLines(
    file,
    warn = FALSE
  )
  
}


#--------------------------------------------------------------
# Extract accession number
#
# Description:
#   Extracts the accession number from a GenBank record.
#
# Arguments:
#   lines : Character vector containing the GenBank record.
#
# Returns:
#   Accession number.
#--------------------------------------------------------------

extract_accession <- function(lines){
  
  if(!is_character_vector(lines)){
    
    stop(
      "'lines' must be a non-empty character vector."
    )
    
  }
  
  accession_line <- grep(
    
    ACCESSION_START,
    
    lines,
    
    value = TRUE
    
  )
  
  if(length(accession_line) != 1){
    
    stop(
      "Unable to identify a unique ACCESSION line."
    )
    
  }
  
  accession <- sub(
    ACCESSION_START,
    "",
    accession_line
  )
  
  return(
    trimws(accession)
  )
  
}

#--------------------------------------------------------------
# Extract organism
#
# Description:
#   Extracts the organism name from a GenBank record.
#
# Arguments:
#   lines : Character vector containing the GenBank record.
#
# Returns:
#   Organism name.
#--------------------------------------------------------------

extract_organism <- function(lines){
  
  if(!is_character_vector(lines)){
    
    stop(
      "'lines' must be a non-empty character vector."
    )
    
  }
  
  organism_line <- grep(
    
    ORGANISM_START,
    
    lines,
    
    value = TRUE
    
  )
  
  if(length(organism_line) != 1){
    
    stop(
      "Unable to identify a unique ORGANISM line."
    )
    
  }
  
  organism <- sub(
    
    ORGANISM_START,
    
    "",
    
    organism_line
    
  )
  
  return(
    trimws(organism)
  )
  
}


#--------------------------------------------------------------
# Extract the FEATURES section
#
# Description:
#   Returns the FEATURES section of a GenBank record exactly as it
#   appears in the source file. Parsing of individual features is
#   handled elsewhere.
#
# Arguments:
#   lines : Character vector containing the GenBank record.
#
# Returns:
#   Character vector containing the FEATURES section.
#--------------------------------------------------------------

extract_features <- function(lines){
  
  if(!is_character_vector(lines)){
    
    stop(
      "'lines' must be a non-empty character vector."
    )
    
  }
  
  extract_feature_section(lines)
  
}


#--------------------------------------------------------------
# Read the genome sequence
#
# Description:
#   Extracts the DNA sequence from the ORIGIN section of a
#   GenBank record and returns it as a single uppercase
#   character string.
#
# Arguments:
#   gb_lines : Character vector containing the GenBank record.
#
# Returns:
#   Character string containing the complete genome sequence.
#--------------------------------------------------------------

read_sequence <- function(gb_lines){
  
  if(!is_character_vector(gb_lines)){
    
    stop(
      "'gb_lines' must be a non-empty character vector."
    )
    
  }
  
  # Locate the ORIGIN line
  origin_start <- grep(
    ORIGIN_START,
    gb_lines
  )
  
  if(length(origin_start) != 1){
    
    stop(
      "Unable to identify a unique ORIGIN section."
    )
    
  }
  
  # Locate the end of the sequence
  origin_end <- grep(
    GENBANK_END,
    gb_lines
  )
  
  if(length(origin_end) != 1){
    
    stop(
      "GenBank record is missing the terminating '//'."
    )
    
  }
  
  # Extract sequence lines
  seq_lines <- gb_lines[(origin_start + 1):(origin_end - 1)]
  
  # Remove position numbers
  seq_lines <- gsub("[[:digit:]]", "", seq_lines)
  
  # Remove spaces
  seq_lines <- gsub("\\s+", "", seq_lines)
  
  # Collapse into one sequence
  sequence <- paste(seq_lines, collapse = "")
  
  # Convert to uppercase
  sequence <- trim_and_upper(sequence)
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Genome sequence contains invalid nucleotide symbols."
    )
    
  }
  
  return(sequence)
  
}