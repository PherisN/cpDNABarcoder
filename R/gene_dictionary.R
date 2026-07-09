###############################################################
##
## gene_dictionary.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Define accepted chloroplast gene names and their
##   annotation metadata.
##
## Responsibilities:
##   - Store accepted gene names
##   - Store accepted aliases
##   - Define preferred GenBank feature types
##   - Define expected molecule types
##   - Define expected copy numbers
##
###############################################################


#==============================================================
# Gene Dictionary
#==============================================================

#--------------------------------------------------------------
# Gene dictionary structure
#
# Description:
#   Each entry defines the accepted representation of a
#   chloroplast gene throughout the package.
#
# Fields:
#   aliases          : Accepted gene names and synonyms.
#   feature_type     : Preferred GenBank annotation.
#   molecule         : Expected biological molecule.
#   expected_copies  : Expected number of genomic copies.
#
# Exported gene dictionary used throughout the package.
#--------------------------------------------------------------

GENE_DICTIONARY <- list(

###############################################################
##
## Protein-coding genes
##
###############################################################  
  PSBA = list(
    
    aliases = c("PSBA"),
    
    feature_type = FEATURE_GENE,
    
    expected_copies = 1
    
  ),
  
  RBCL = list(
    
    aliases = c("RBCL"),
    
    feature_type = FEATURE_GENE,
    
    expected_copies = 1
    
  ),
  
  MATK = list(
    
    aliases = c("MATK"),
    
    feature_type = FEATURE_GENE,
    
    expected_copies = 1
    
  ),
  
###############################################################
##
## Transfer RNA genes
##
###############################################################  
  
  TRNH = list(
    
    aliases = c(
      "TRNH",
      "TRNH-GUG",
      "TRNH-GTG"
    ),
    
    feature_type = FEATURE_GENE,
    
    expected_copies = 1
    
  )

###############################################################
##
## Ribosomal RNA genes
##
###############################################################

# Reserved for future additions.


)


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate a gene dictionary
#
# Description:
#   Checks that every gene dictionary entry contains the
#   required fields and that each field contains a valid
#   value.
#
# Arguments:
#   dictionary : Gene dictionary.
#
# Returns:
#   Invisibly returns the validated dictionary.
#--------------------------------------------------------------

validate_gene_dictionary <- function(dictionary){
  
  if(!is.list(dictionary)){
    
    stop(
      "'dictionary' must be a list."
    )
    
  }
  
  required <- c(
    
    "aliases",
    
    "feature_type",
    
    "expected_copies"
    
  )
  
  valid_features <- c(
    
    FEATURE_GENE,
    
    FEATURE_CDS,
    
    FEATURE_TRNA,
    
    FEATURE_RRNA
    
  )
  
  for(name in names(dictionary)){
    
    entry <- dictionary[[name]]
    
    #-------------------------------------------
    # Required fields
    #-------------------------------------------
    
    missing <- required[
      !required %in% names(entry)
    ]
    
    if(length(missing) > 0){
      
      stop(
        
        "Gene '",
        
        name,
        
        "' is missing field(s): ",
        
        paste(
          missing,
          collapse = ", "
        )
        
      )
      
    }
    
    #-------------------------------------------
    # Aliases
    #-------------------------------------------
    
    if(
      !is.character(entry$aliases) ||
      length(entry$aliases) == 0 ||
      any(is.na(entry$aliases))
    ){
      
      stop(
        
        "Gene '",
        
        name,
        
        "' must define one or more aliases."
        
      )
      
    }
    
    #-------------------------------------------
    # Feature type
    #-------------------------------------------
    
    if(
      !is_single_string(entry$feature_type) ||
      !(entry$feature_type %in% valid_features)
    ){
      
      stop(
        
        "Gene '",
        
        name,
        
        "' has an invalid feature_type."
        
      )
      
    }
    
    #-------------------------------------------
    # Expected copy number
    #-------------------------------------------
    
    if(
      !is.numeric(entry$expected_copies) ||
      length(entry$expected_copies) != 1 ||
      is.na(entry$expected_copies) ||
      entry$expected_copies < 1 ||
      entry$expected_copies != as.integer(entry$expected_copies)
    ){
      
      stop(
        
        "Gene '",
        
        name,
        
        "' has an invalid expected_copies value."
        
      )
      
    }
    
  }
  
  invisible(dictionary)
  
}

#--------------------------------------------------------------
# Validate the package gene dictionary
#--------------------------------------------------------------

validate_gene_dictionary(
  GENE_DICTIONARY
)
