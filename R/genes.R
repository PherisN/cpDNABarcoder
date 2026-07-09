###############################################################
##
## genes.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Locate and validate annotated genes in chloroplast
##   genomes.
##
## Responsibilities:
##   - Standardize gene names
##   - Locate annotated genes
##   - Select preferred GenBank features
##   - Validate selected annotations
##
###############################################################

#==============================================================
# Gene Utilities
#==============================================================

#--------------------------------------------------------------
# Standardize a gene name
#
# Description:
#   Converts a gene name to the standard representation used
#   throughout the package. Gene names are converted to
#   uppercase, leading and trailing whitespace is removed, and
#   underscores are replaced with hyphens.
#
# Arguments:
#   gene : Gene name.
#
# Returns:
#   Standardized gene name.
#--------------------------------------------------------------

standardize_gene_name <- function(gene){
  
  if(!is_single_string(gene)){
    
    stop(
      "'gene' must be a single character string."
    )
    
  }
  
  gene <- trim_and_upper(gene)
  
  gene <- gsub(
    "_",
    "-",
    gene
  )
  
  return(gene)
  
}


#==============================================================
# Gene Lookup
#==============================================================

#--------------------------------------------------------------
# Find annotated gene features
#
# Description:
#   Searches the parsed GenBank feature table using all
#   accepted aliases defined in the package gene dictionary.
#
# Arguments:
#   feature_table   : Parsed feature table.
#   gene_name       : Gene to locate.
#   gene_dictionary : Gene dictionary used for matching.
#
# Returns:
#   A data frame containing all matching feature records.
#--------------------------------------------------------------

find_gene <- function(feature_table,
                      gene_name,
                      gene_dictionary = GENE_DICTIONARY){
  
  if(!is.data.frame(feature_table)){
    
    stop(
      "'feature_table' must be a data.frame."
    )
    
  }
  
  required <- c(
    "gene",
    "type",
    "start",
    "end",
    "strand"
  )
  
  missing <- setdiff(
    required,
    names(feature_table)
  )
  
  if(length(missing) > 0){
    
    stop(
      "Feature table is missing required columns: ",
      paste(
        missing,
        collapse = ", "
      )
    )
    
  }
  
  gene_name <- standardize_gene_name(gene_name)
  
  if(!gene_name %in% names(gene_dictionary)){
    
    stop(
      "Gene '",
      gene_name,
      "' is not defined in GENE_DICTIONARY."
    )
    
  }
  
  gene_info <- gene_dictionary[[gene_name]]
  
  aliases <- vapply(
    gene_info$aliases,
    standardize_gene_name,
    character(1)
  )
  
  genes <- rep(NA_character_, nrow(feature_table))
  
  has_gene <- !is.na(feature_table$gene)
  
  genes[has_gene] <- vapply(
    feature_table$gene[has_gene],
    standardize_gene_name,
    character(1)
  )
  
  matches <- feature_table[
    genes %in% aliases,
    ,
    drop = FALSE
  ]
  
  if(nrow(matches) == 0){
    
    warning(
      paste0(
        "No annotations found for ",
        gene_name
      )
    )
    
  }
  
  return(matches)
  
}


#--------------------------------------------------------------
# Select the preferred gene feature
#
# Description:
#   Selects the preferred GenBank feature for a gene according
#   to the feature type specified in the gene dictionary.
#
# Arguments:
#   matches         : Matching feature records.
#   gene_name       : Target gene.
#   gene_dictionary : Gene dictionary.
#
# Returns:
#   The preferred feature record(s).
#--------------------------------------------------------------

select_gene_feature <- function(matches,
                                gene_name,
                                gene_dictionary = GENE_DICTIONARY){
  
  gene_name <- standardize_gene_name(gene_name)
  
  if(!gene_name %in% names(gene_dictionary)){
    
    stop(
      "Unknown gene: ",
      gene_name
    )
    
  }
  
  gene_info <- gene_dictionary[[gene_name]]
  
  preferred <- gene_info$feature_type
  
  feature <- matches[
    matches$type == preferred,
  ]
  
  if(nrow(feature) == 0){
    
    warning(
      paste(
        "Preferred feature",
        preferred,
        "not found."
      )
    )
    
    return(matches)
    
  }
  
  if(nrow(feature) > 1){
    
    warning(
      paste(
        "Multiple",
        preferred,
        "features found."
      )
    )
    
  }
  
  return(feature)
  
}


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate a selected genomic feature
#
# Description:
#   Verifies that a selected feature contains usable genomic
#   coordinates and strand information. If multiple feature
#   records are supplied, the first record is retained after
#   issuing a warning.
#
# Arguments:
#   feature         : Selected feature record(s).
#   gene_name       : Target gene.
#   gene_dictionary : Gene dictionary.
#
# Returns:
#   A validated feature record.
#--------------------------------------------------------------

validate_feature <- function(feature,
                             gene_name){
  
  # Empty result?
  if (nrow(feature) == 0) {
    
    stop(
      paste0(
        "No annotation found for gene: ",
        gene_name
      )
    )
    
  }
  
  # Multiple rows?
  if (nrow(feature) > 1) {
    
    warning(
      paste(
        "Multiple annotations found for",
        gene_name,
        "; using the first record."
      )
    )
    
    feature <- feature[1, ]
    
  }
  
  # Coordinates available?
  if(
    !is_valid_coordinate(feature$start) ||
    !is_valid_coordinate(feature$end)
  ){
    
    stop(
      paste0(
        "Coordinates could not be parsed for ",
        gene_name
      )
    )
    
  }
  
  # Strand known?
  if (!(feature$strand %in% c(
    STRAND_FORWARD,
    STRAND_REVERSE
  ))) {
    
    stop(
      paste0(
        "Unknown strand for ",
        gene_name
      )
    )
    
  }
  
  if(
    is.na(feature$type) ||
    !nzchar(feature$type)
  ){
    
    stop(
      "Feature type is missing."
    )
    
  }
  
  return(feature)
  
}

#==============================================================
# High-level Interface
#==============================================================

#--------------------------------------------------------------
# Retrieve a validated gene feature
#
# Description:
#   Retrieves a validated genomic feature corresponding to the
#   requested gene by locating candidate annotations,
#   selecting the preferred feature, and validating the final
#   result.
#
# Arguments:
#   feature_table   : Parsed GenBank feature table.
#   gene_name       : Target gene.
#   gene_dictionary : Gene dictionary.
#
# Returns:
#   A validated genomic feature.
#--------------------------------------------------------------

get_gene_feature <- function(feature_table,
                             gene_name,
                             gene_dictionary = GENE_DICTIONARY){
  
  matches <- find_gene(
    feature_table,
    gene_name,
    gene_dictionary
  )
  
  selected <- select_gene_feature(
    matches,
    gene_name,
    gene_dictionary
  )
  
  validated <- validate_feature(
    selected,
    gene_name
  )
  
  return(validated)
  
}
