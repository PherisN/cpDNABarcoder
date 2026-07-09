###############################################################
##
## parser.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Parse the FEATURES section of GenBank records into
##   structured R objects for downstream marker extraction.
##
## Responsibilities:
##   - Extract the FEATURES section
##   - Identify feature boundaries
##   - Parse feature locations
##   - Parse feature qualifiers
##   - Construct a feature table
##
###############################################################

#==============================================================
# Feature Utilities
#==============================================================

#--------------------------------------------------------------
# Identify the start of a GenBank feature
#
# Description:
#   Determines whether a line begins a new GenBank feature.
#
# Arguments:
#   line : Character string.
#
# Returns:
#   TRUE if the line starts a feature, FALSE otherwise.
#--------------------------------------------------------------

is_feature_start <- function(line){
  
  if(!is_single_string(line)){
    
    stop(
      "'line' must be a single character string."
    )
    
  }
  
  grepl("^     [A-Za-z]", line)
  
}

#--------------------------------------------------------------
# Extract the FEATURES section
#
# Description:
#   Extracts the FEATURES section from a GenBank record.
#
# Arguments:
#   gb_lines : Character vector containing the GenBank record.
#
# Returns:
#   Character vector containing only the FEATURES section.
#--------------------------------------------------------------

extract_feature_section <- function(gb_lines){
  
  if(!is_character_vector(gb_lines)){
    
    stop(
      "'gb_lines' must be a non-empty character vector."
    )
    
  }
  
  start <- grep(
    GENBANK_FEATURE_START,
    gb_lines
  )
  
  if(length(start) != 1){
    
    stop(
      "Unable to identify a unique FEATURES section."
    )
    
  }
  
  end <- grep(
    ORIGIN_START,
    gb_lines
  )
  
  if(length(end) != 1){
    
    stop(
      "Unable to identify a unique ORIGIN section."
    )
    
  }
  
  features <- gb_lines[(start + 1):(end - 1)]
  
  return(features)
  
}


#==============================================================
# Feature Constructors
#==============================================================

#--------------------------------------------------------------
# Create an empty feature record
#
# Description:
#   Creates an empty feature object with the fields required
#   during GenBank feature parsing.
#
# Returns:
#   A feature record represented as a named list.
#--------------------------------------------------------------

new_feature <- function() {
  
  list(
    type = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    strand = NA_character_,
    location_type = NA_character_,
    location = NA_character_,     
    gene = NA_character_,
    product = NA_character_,
    locus_tag = NA_character_,
    protein_id = NA_character_,
    note = NA_character_
  )
  
}

#--------------------------------------------------------------
# Parse a GenBank feature header
#
# Description:
#   Parses the first line of a GenBank feature block into
#   feature type and genomic location.
#
# Arguments:
#   line : Character string containing the feature header.
#
# Returns:
#   A list containing:
#     - type
#     - location
#--------------------------------------------------------------

parse_feature_header <- function(line){
  
  if(!is_single_string(line)){
    
    stop(
      "'line' must be a single character string."
    )
    
  }
  
  line <- trimws(line)
  
  parts <- strsplit(line, "\\s+")[[1]]
  
  if(length(parts) < 2){
    
    stop(
      "Invalid feature header."
    )
    
  }
  
  type <- parts[1]
  
  location <- paste(parts[-1], collapse = " ")
  
  list(
    
    type = type,
    
    location = location
    
  )
  
}

#--------------------------------------------------------------
# Parse a GenBank location string
#
# Description:
#   Parses GenBank genomic location strings.
#
# Supported location types:
#   - simple
#   - complement(...)
#   - join(...)
#   - order(...)
#
# Notes:
#   Joined and ordered locations currently return missing
#   coordinates (NA) while preserving the location type.
#
# Arguments:
#   location : GenBank location string.
#
# Returns:
#   A parsed location record.
#--------------------------------------------------------------

parse_location <- function(location){
  
  if(!is_single_string(location)){
    
    stop(
      "'location' must be a single character string."
    )
    
  }
  
  location <- trimws(location)
  
  strand <- STRAND_FORWARD
  
  if(grepl("^complement\\(", location)){
    
    strand <- STRAND_REVERSE
    
    location <- sub("^complement\\(", "", location)
    
    location <- sub("\\)$", "", location)
    
  }
  
  location_type <- LOCATION_SIMPLE
  
  if(grepl("^join\\(", location)){
    
    return(list(
      
      start = NA_integer_,
      
      end = NA_integer_,
      
      strand = strand,
      
      location_type = LOCATION_JOIN
      
    ))
    
  }
  
  if(grepl("^order\\(", location)){
    
    return(list(
      
      start = NA_integer_,
      
      end = NA_integer_,
      
      strand = strand,
      
      location_type = LOCATION_ORDER
      
    ))
    
  }
  
  location <- gsub("[<>]", "", location)
  
  coords <- strsplit(
    location,
    "\\.\\."
  )[[1]]
  
  if(length(coords) != 2){
    
    stop(
      "Invalid GenBank location string."
    )
    
  }
  
  start <- as.integer(coords[1])
  
  end <- as.integer(coords[2])
  
  if(is.na(start) || is.na(end)){
    
    stop(
      "Unable to parse genomic coordinates."
    )
    
  }
  
  list(
    
    start = start,
    
    end = end,
    
    strand = strand,
    
    location_type = location_type,
    
    raw_location = location
    
  )
  
}

#--------------------------------------------------------------
# Parse a GenBank feature qualifier
#
# Description:
#   Parses a qualifier line from a GenBank feature block.
#
# Arguments:
#   line : Character string representing one qualifier.
#
# Returns:
#   A list containing:
#     - key   : qualifier name
#     - value : qualifier value
#--------------------------------------------------------------

parse_qualifier <- function(line){
  
  if(!is_single_string(line)){
    
    stop(
      "'line' must be a single character string."
    )
    
  }
  
  # Remove leading whitespace
  line <- trimws(line)
  
  # Remove leading "/"
  line <- sub("^/", "", line)
  
  # Qualifiers without an assigned value
  if(!grepl("=", line)){
    
    return(
      list(
        key = line,
        value = TRUE
      )
    )
    
  }
  
  # Split key and value
  parts <- strsplit(
    line,
    "=",
    fixed = TRUE
  )[[1]]
  
  key <- parts[1]
  
  value <- paste(parts[-1], collapse = "=")
  
  # Remove quotation marks
  value <- gsub('"', "", value)
  
  return(
    list(
      key = key,
      value = value
    )
  )
  
}

#--------------------------------------------------------------
# Split feature blocks
#
# Description:
#   Splits the FEATURES section into individual GenBank feature
#   blocks.
#
# Arguments:
#   features : Character vector containing the FEATURES
#              section.
#
# Returns:
#   A list of feature blocks.
#--------------------------------------------------------------

split_feature_blocks <- function(features){
  
  if(!is_character_vector(features)){
    
    stop(
      "'features' must be a non-empty character vector."
    )
    
  }
  
  blocks <- list()
  
  current_block <- character()
  
  for(line in features){
    
    if(is_feature_start(line)){
      
      if(length(current_block) > 0){
        
        blocks[[length(blocks) + 1]] <- current_block
        
      }
      
      current_block <- line
      
    } else{
      
      current_block <- c(current_block, line)
      
    }
    
  }
  
  # Add the last feature
  if(length(current_block) > 0){
    
    blocks[[length(blocks)+1]] <- current_block
    
  }
  
  if(length(blocks) == 0){
    
    warning(
      "No feature blocks were identified."
    )
    
  }
  
  return(blocks)
  
}

#--------------------------------------------------------------
# Parse a single GenBank feature block
#
# Description:
#   Parses one feature block into a structured feature record.
#
# Arguments:
#   block : Character vector containing one feature block.
#
# Returns:
#   A parsed feature (list).
#--------------------------------------------------------------

parse_feature <- function(block){
  
  if(!is_character_vector(block)){
    
    stop(
      "'block' must be a non-empty character vector."
    )
    
  }
  
  feature <- new_feature()
  
  # Parse feature header
  header <- parse_feature_header(block[1])
  
  feature$type <- header$type
  feature$location <- header$location
  
  # Parse genomic location
  loc <- parse_location(header$location)
  
  feature$start <- loc$start
  feature$end <- loc$end
  feature$strand <- loc$strand
  feature$location_type <- loc$location_type
  
  # Parse qualifiers
  if(length(block) >= 2){
    
    for(i in seq.int(2, length(block))){
      
      qualifier <- parse_qualifier(block[i])
      
      if(qualifier$key %in% names(feature)){
        
        feature[[qualifier$key]] <- qualifier$value
        
      }
      
    }
    
  }
  
  return(feature)
  
}

#==============================================================
# High-level Parser
#==============================================================

#--------------------------------------------------------------
# Parse all FEATURES from a GenBank record
#
# Description:
#   Extracts the FEATURES section from a GenBank file,
#   splits it into individual feature blocks,
#   parses each block, and returns a feature table.
#
# Arguments:
#   gb_lines : Character vector containing the GenBank file.
#
# Returns:
#   A data.frame where each row represents one parsed feature.
#--------------------------------------------------------------

parse_features <- function(gb_lines){
  
  if(!is_character_vector(gb_lines)){
    
    stop(
      "'gb_lines' must be a non-empty character vector."
    )
    
  }
  
  # Extract the FEATURES section
  features <- extract_feature_section(gb_lines)
  
  # Split into individual feature blocks
  blocks <- split_feature_blocks(features)
  
  if(length(blocks) == 0){
    
    stop(
      "No feature blocks found."
      )
    
  }
  
  # Parse each feature block
  feature_list <- lapply(
    blocks,
    parse_feature
  )
  
  # Combine into a single table
  feature_table <- do.call(
    rbind,
    lapply(
      feature_list,
      as.data.frame
    )
  )
  
  # Clean row names
  rownames(feature_table) <- NULL
  
  return(feature_table)
  
}