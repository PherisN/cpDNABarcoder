###############################################################
##
## marker_dictionary.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Define supported DNA barcode markers and the metadata
##   required to extract them.
##
## Responsibilities:
##   - Define supported barcode markers
##   - Classify marker types
##   - Store extraction metadata
##   - Support marker dispatch
##
###############################################################

#==============================================================
# Marker Dictionary
#==============================================================

#--------------------------------------------------------------
# Marker dictionary structure
#
# Description:

#   Each entry defines the information required to identify
#   and extract a supported DNA barcode marker.
#
# Fields:
#   class        : Marker class.
#   marker       : Standard marker name.
#   description  : Human-readable description.
#   type         : Genome type.
#   gene         : Target gene.
#   left_gene    : Left flanking gene.
#   right_gene   : Right flanking gene.
#   locus        : Nuclear locus.
#--------------------------------------------------------------

MARKER_DICTIONARY <- list(
  ###############################################################
  ##
  ## Intergenic markers
  ##
  ###############################################################
  
  PSBA_TRNH = list(
    
    class = MARKER_CLASS_INTERGENIC,
    
    marker = "psbA-trnH",
    
    description =
      "Intergenic spacer between psbA and trnH",
    
    type = "plastid",
    
    left_gene = "PSBA",
    
    right_gene = "TRNH"
    
  ),
  
  
  ###############################################################
  ##
  ## Coding gene markers
  ##
  ###############################################################
  
  RBCL = list(
    
    class = MARKER_CLASS_GENE,
    
    marker = "rbcL",
    
    description =
      "Large subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase (RuBisCO)",
    
    type = "plastid",
    
    gene = "RBCL"
    
  ),
  
  
  MATK = list(
    
    class = MARKER_CLASS_GENE,
    
    marker = "matK",
    
    description =
      "Maturase K",
    
    type = "plastid",
    
    gene = "MATK"
    
  )
  
)
  ###############################################################
  ##
  ## Introns markers
  ##
  ###############################################################
  
  # Reserved for future markers.
  
  ###############################################################
  ##
  ## Mitochondrial markers
  ##
  ###############################################################
  
  # Reserved for future markers.


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate a marker dictionary
#
# Description:
#   Checks that every marker definition contains the required
#   fields for its marker class and that each field contains
#   a valid value.
#
# Arguments:
#   dictionary : Marker dictionary.
#
# Returns:
#   Invisibly returns the validated dictionary.
#--------------------------------------------------------------

validate_marker_dictionary <- function(dictionary) {
  
  if(!is.list(dictionary)) {
    stop("'dictionary' must be a list.")
    
  }
  
  valid_classes <- c(MARKER_CLASS_GENE,
                     
                     MARKER_CLASS_INTERGENIC,
                     
                     MARKER_CLASS_NUCLEAR)
  
  valid_types <- c("plastid", "nuclear")
  
  for (name in names(dictionary)) {
    marker <- dictionary[[name]]
    
    #-------------------------------------------
    # Required common fields
    #-------------------------------------------
    
    required <- c("class", "marker", "description", "type")
    
    missing <- required[!required %in% names(marker)]
    
    if(length(missing) > 0) {
      stop("Marker '",
           
           name,
           
           "' is missing field(s): ",
           
           paste(missing, collapse = ", "))
      
    }
    
    #-------------------------------------------
    # Marker class
    #-------------------------------------------
    
    if(!is_single_string(marker$class) ||
        !(marker$class %in% valid_classes)) {
      stop("Marker '", name, "' has an invalid marker class.")
      
    }
    
    #-------------------------------------------
    # Marker name
    #-------------------------------------------
    
    if(!is_single_string(marker$marker)) {
      stop("Marker '", name, "' has an invalid marker name.")
      
    }
    
    #-------------------------------------------
    # Description
    #-------------------------------------------
    
    if(!is_single_string(marker$description)) {
      stop("Marker '", name, "' has an invalid description.")
      
    }
    
    #-------------------------------------------
    # Genome type
    #-------------------------------------------
    
    if(!is_single_string(marker$type) ||
        !(marker$type %in% valid_types)) {
      stop("Marker '", name, "' has an invalid genome type.")
      
    }
    
    #-------------------------------------------
    # Class-specific validation
    #-------------------------------------------
    
    if(marker$class == MARKER_CLASS_GENE) {
      if(!is_single_string(marker$gene)) {
        stop("Marker '", name, "' has an invalid gene.")
        
      }
      
    }
    
    if(marker$class == MARKER_CLASS_INTERGENIC) {
      if(!is_single_string(marker$left_gene)) {
        stop("Marker '", name, "' has an invalid left_gene.")
        
      }
      
      if(!is_single_string(marker$right_gene)) {
        stop("Marker '", name, "' has an invalid right_gene.")
        
      }
      
    }
    
    if(marker$class == MARKER_CLASS_NUCLEAR) {
      if(!is_single_string(marker$locus)) {
        stop("Marker '", name, "' has an invalid locus.")
        
      }
      
    }
    
  }
  
  invisible(dictionary)
  
}

#--------------------------------------------------------------
# Validate the package marker dictionary
#--------------------------------------------------------------

validate_marker_dictionary(MARKER_DICTIONARY)