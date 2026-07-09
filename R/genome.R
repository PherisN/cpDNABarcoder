###############################################################
##
## genome.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Construct and manage cpGenome objects.
##
## Responsibilities:
##   - Read GenBank genomes
##   - Store genome metadata
##   - Store genome sequence
##   - Store parsed feature annotations
##
###############################################################

#==============================================================
# Genome Constructors
#==============================================================

#--------------------------------------------------------------
# Create a cpGenome object
#
# Description:
#   Constructs a chloroplast genome object containing genome
#   metadata, DNA sequence, and parsed feature annotations.
#
# Arguments:
#   metadata      : Genome metadata list.
#   sequence      : Character string containing the genome
#                   sequence.
#   feature_table : Parsed GenBank feature table.
#
# Returns:
#   A validated cpGenome object.
#--------------------------------------------------------------

new_cpGenome <- function(metadata,
                         sequence,
                         feature_table){
  
  if(!is.list(metadata)){
    
    stop(
      "'metadata' must be a list."
    )
    
  }
  
  if(!is_single_string(sequence)){
    
    stop(
      "'sequence' must be a single character string."
    )
    
  }
  
  if(!is.data.frame(feature_table)){
    
    stop(
      "'feature_table' must be a data.frame."
    )
    
  }
  
  genome <- list(
    
    metadata = metadata,
    
    sequence = sequence,
    
    feature_table = feature_table
    
  )
  
  class(genome) <- "cpGenome"
  
  return(
    validate_cpGenome(genome)
  )
  
}


#--------------------------------------------------------------
# Validate a cpGenome object
#
# Description:
#   Checks that a cpGenome object contains all required
#   components.
#
# Arguments:
#   genome : A cpGenome object.
#
# Returns:
#   The validated cpGenome object.
#--------------------------------------------------------------

validate_cpGenome <- function(genome){
  
  if(!inherits(genome, "cpGenome")){
    
    stop(
      "Object is not a cpGenome."
    )
    
  }
  
  if(!is.list(genome$metadata)){
    
    stop(
      "Genome metadata must be a list."
    )
    
  }
  
  if(!is.character(genome$sequence) ||
     length(genome$sequence) != 1){
    
    stop(
      "Genome sequence must be a single character string."
    )
    
  }
  
  if(nchar(genome$sequence) == 0){
    
    stop(
      "Genome sequence is empty."
    )
    
  }
  
  if(!is_valid_dna(genome$sequence)){
    
    stop(
      "Genome sequence contains invalid DNA symbols."
    )
    
  }
  
  if(!is.data.frame(genome$feature_table)){
    
    stop(
      "Feature table must be a data.frame."
    )
    
  }
  
  return(genome)
  
}


#--------------------------------------------------------------
# Assert that an object is a cpGenome
#
# Description:
#   Stops execution if the supplied object is not a cpGenome.
#
# Arguments:
#   genome : Object to test.
#
# Returns:
#   Invisibly returns TRUE.
#--------------------------------------------------------------

assert_cpGenome <- function(genome){
  
  if(!inherits(genome, "cpGenome")){
    
    stop(
      "Object is not a cpGenome."
    )
    
  }
  
  invisible(genome)
  
}


#==============================================================
# Metadata
#==============================================================

#--------------------------------------------------------------
# Read genome metadata
#
# Description:
#   Extracts genome metadata from a GenBank record.
#
# Arguments:
#   gb_lines : Character vector containing a GenBank record.
#
# Returns:
#   A metadata list.
#--------------------------------------------------------------

read_metadata <- function(gb_lines){
  
  if(!is_character_vector(gb_lines)){
    
    stop(
      "'gb_lines' must be a non-empty character vector."
    )
    
  }
  
  accession_line <- grep(
    ACCESSION_START,
    gb_lines,
    value = TRUE
  )
  
  if(length(accession_line) != 1){
    
    stop(
      "Unable to identify a unique ACCESSION line."
    )
    
  }
  
  accession <- trimws(
    sub(
      ACCESSION_START,
      "",
      accession_line
    )
  )
  
  version_line <- grep(
    VERSION_START,
    gb_lines,
    value = TRUE
  )
  
  if(length(version_line) != 1){
    
    stop(
      "Unable to identify a unique VERSION line."
    )
    
  }
  
  version <- trimws(
    sub(
      VERSION_START,
      "",
      version_line
    )
  )
  
  definition_line <- grep(
    DEFINITION_START,
    gb_lines,
    value = TRUE
  )
  
  if(length(definition_line) != 1){
    
    stop(
      "Unable to identify a unique DEFINITION line."
    )
    
  }
  
  definition <- trimws(
    sub(
      DEFINITION_START,
      "",
      definition_line
    )
  )
  
  organism_line <- grep(
    ORGANISM_START,
    gb_lines,
    value = TRUE
  )
  
  if(length(organism_line) != 1){
    
    stop(
      "Unable to identify a unique ORGANISM line."
    )
    
  }
  
  organism <- trimws(
    sub(
      ORGANISM_START,
      "",
      organism_line
    )
  )
  
  locus_line <- grep(
    LOCUS_START,
    gb_lines,
    value = TRUE
  )
  
  if(length(locus_line) != 1){
    
    stop(
      "Unable to identify a unique LOCUS line."
    )
    
  }
  
  genome_length <- as.integer(
    strsplit(
      trimws(locus_line),
      "\\s+"
    )[[1]][3]
  )
  
  if(is.na(genome_length)){
    
    stop(
      "Unable to determine genome length."
    )
    
  }
  
  list(
    
    accession = accession,
    
    version = version,
    
    organism = organism,
    
    definition = definition,
    
    genome_length = genome_length
    
  )
  
}


#==============================================================
# High-level Readers
#==============================================================

#--------------------------------------------------------------
# Read a chloroplast genome
#
# Description:
#   Reads a GenBank file and constructs a validated cpGenome
#   object.
#
# Arguments:
#   file : Path to a GenBank file.
#
# Returns:
#   A validated cpGenome object.
#--------------------------------------------------------------

read_cpGenome <- function(file){
    
    gb <- read_genbank(file)
    
    metadata <- read_metadata(gb)
    
    sequence <- read_sequence(gb)
    
    feature_table <- parse_features(gb)
    
    genome <- new_cpGenome(
      
      metadata = metadata,
      
      sequence = sequence,
      
      feature_table = feature_table
      
    )
    
    return(
      validate_cpGenome(genome)
    )
    
}


#==============================================================
# Accessors
#==============================================================

#--------------------------------------------------------------
# cpGenome accessors
#--------------------------------------------------------------

get_metadata.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata)
  
}

get_sequence.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$sequence)
  
}

get_feature_table.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$feature_table)
  
}

get_accession.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata$accession)
  
}

get_version.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata$version)
  
}

get_definition.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata$definition)
  
}

get_organism.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata$organism)
  
}

get_genome_length.cpGenome <- function(x){
  
  assert_cpGenome(x)
  
  return(x$metadata$genome_length)
  
}


#==============================================================
# Methods
#==============================================================

#--------------------------------------------------------------
# Print a cpGenome object
#
# Description:
#   Displays a concise summary of a cpGenome object.
#
# Arguments:
#   x : A cpGenome object.
#
# Returns:
#   Invisibly returns the input object.
#--------------------------------------------------------------

print.cpGenome <- function(x, ...){
  
  assert_cpGenome(x)
  
  cat("\n")
  
  cat("cpGenome Object\n")
  
  cat("---------------------------\n")
  
  cat(
    "Organism :",
    x$metadata$organism,
    "\n"
  )
  
  cat(
    "Accession:",
    x$metadata$accession,
    "\n"
  )
  
  cat(
    "Genome   :",
    x$metadata$genome_length,
    "bp\n"
  )
  
  cat(
    "Sequence :",
    nchar(x$sequence),
    "bp\n"
  )
  
  cat(
    "Features :",
    nrow(x$feature_table),
    "annotations\n"
  )
  
  invisible(x)
  
}