###############################################################
##
## pipeline_options.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Define and validate configurable options controlling
##   DNA barcode extraction workflows.
##
## Responsibilities:
##   - Define default pipeline options
##   - Validate pipeline option values
##   - Construct customized option sets
##
###############################################################

#==============================================================
# Pipeline Options
#==============================================================

#--------------------------------------------------------------
# Pipeline options
#
# Description:
#   Default configuration controlling DNA barcode extraction
#   and export workflows.
#
# Fields:
#   reverse_complement : Return reverse complement sequences.
#   include_flanks     : Include flanking sequence.
#   flank_length       : Number of flanking bases.
#   allow_wraparound   : Permit extraction across the genome
#                        origin.
#   export_coordinates : Include genomic coordinates in output.
#   verbose            : Print informative messages.
#--------------------------------------------------------------

PIPELINE_OPTIONS <- list(
  
  #############################################################
  ##
  ## Sequence extraction
  ##
  #############################################################
  
  reverse_complement = FALSE,
  
  include_flanks = FALSE,
  
  flank_length = 0,
  
  allow_wraparound = TRUE,
  
  #############################################################
  ##
  ## Output
  ##
  #############################################################
  
  export_coordinates = TRUE,
  
  #############################################################
  ##
  ## Logging
  ##
  #############################################################
  
  verbose = TRUE
  
)


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate pipeline options
#
# Description:
#   Checks that a pipeline options object contains all required
#   settings and that each option has the expected data type.
#
# Arguments:
#   options : Pipeline options list.
#
# Returns:
#   The validated pipeline options list.
#--------------------------------------------------------------

validate_pipeline_options <- function(options){
  
  #---------------------------------------------
  # Validate object
  #---------------------------------------------
  
  if(!is.list(options)){
    
    stop(
      "Pipeline options must be supplied as a list."
    )
    
  }
  
  #---------------------------------------------
  # Validate required fields
  #---------------------------------------------
  
  required <- names(PIPELINE_OPTIONS)
  
  missing <- setdiff(
    required,
    names(options)
  )
  
  if(length(missing) > 0){
    
    stop(
      
      "Missing pipeline option(s): ",
      
      paste(
        missing,
        collapse = ", "
      )
      
    )
    
  }
  
  #---------------------------------------------
  # Reject unknown fields
  #---------------------------------------------
  
  unknown <- setdiff(
    names(options),
    names(PIPELINE_OPTIONS)
  )
  
  if(length(unknown) > 0){
    
    stop(
      
      "Unknown pipeline option(s): ",
      
      paste(
        unknown,
        collapse = ", "
      )
      
    )
    
  }
  
  #---------------------------------------------
  # Validate option types
  #---------------------------------------------
  
  if(
    !is.logical(options$reverse_complement) ||
    length(options$reverse_complement) != 1
  ){
    
    stop(
      "'reverse_complement' must be TRUE or FALSE."
    )
    
  }
  
  if(
    !is.logical(options$include_flanks) ||
    length(options$include_flanks) != 1
  ){
    
    stop(
      "'include_flanks' must be TRUE or FALSE."
    )
    
  }
  
  if(
    !is.numeric(options$flank_length) ||
    length(options$flank_length) != 1 ||
    is.na(options$flank_length) ||
    options$flank_length < 0 ||
    options$flank_length != as.integer(options$flank_length)
  ){
    
    stop(
      "'flank_length' must be a non-negative integer."
    )
    
  }
  
  if(
    !is.logical(options$allow_wraparound) ||
    length(options$allow_wraparound) != 1
  ){
    
    stop(
      "'allow_wraparound' must be TRUE or FALSE."
    )
    
  }
  
  if(
    !is.logical(options$export_coordinates) ||
    length(options$export_coordinates) != 1
  ){
    
    stop(
      "'export_coordinates' must be TRUE or FALSE."
    )
    
  }
  
  if(
    !is.logical(options$verbose) ||
    length(options$verbose) != 1
  ){
    
    stop(
      "'verbose' must be TRUE or FALSE."
    )
    
  }
  
  return(options)
  
}


#==============================================================
# Constructors
#==============================================================

#--------------------------------------------------------------
# Create pipeline options
#
# Description:
#   Creates a validated pipeline options object by overriding 
#   one or more default settings while preserving unspecified 
#   defaults.
#
# Arguments:
#   ... : Named option values to override.
#
# Returns:
#   A validated pipeline options list.
#--------------------------------------------------------------

new_pipeline_options <- function(...){
  
  options <- modifyList(
    
    PIPELINE_OPTIONS,
    
    list(...)
    
  )
  
  return(
    
    validate_pipeline_options(options)
    
  )
  
}


#--------------------------------------------------------------
# Validate package default pipeline options
#--------------------------------------------------------------

validate_pipeline_options(
  
  PIPELINE_OPTIONS
  
)