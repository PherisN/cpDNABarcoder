###############################################################
##
## interval.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Construct and manage genomic interval objects used for
##   intergenic DNA barcode extraction.
##
## Responsibilities:
##   - Construct Interval objects
##   - Validate Interval objects
##   - Provide Interval accessors
##   - Display Interval summaries
##
###############################################################

#==============================================================
# Constructors
#==============================================================

#--------------------------------------------------------------
# Create an Interval object
#
# Description:
#   Creates an S3 object representing the genomic interval
#   between two annotated features.
#
# Arguments:
#   start         : Start coordinate.
#   end           : End coordinate.
#   wraps_origin  : Logical indicating whether the interval
#                   crosses the origin of a circular genome.
#   left_gene     : Name of the left flanking gene.
#   right_gene    : Name of the right flanking gene.
#
# Returns:
#   A validated Interval object.
#--------------------------------------------------------------

new_interval <- function(start,
                         end,
                         wraps_origin,
                         left_gene,
                         right_gene){
  
  if(!is_valid_coordinate(start)){
    stop("'start' must be a valid genomic coordinate.")
  }
  
  if(!is_valid_coordinate(end)){
    stop("'end' must be a valid genomic coordinate.")
  }
  
  if(!is.logical(wraps_origin) ||
     length(wraps_origin) != 1){
    stop("'wraps_origin' must be TRUE or FALSE.")
  }
  
  if(!is_single_string(left_gene)){
    stop("'left_gene' must be a single character string.")
  }
  
  if(!is_single_string(right_gene)){
    stop("'right_gene' must be a single character string.")
  }
  
  width <- if(wraps_origin){
    
    NA_integer_
    
  }else{
    
    end - start + 1
    
  }
  
  interval <- list(
    
    start = start,
    
    end = end,
    
    width = width,
    
    wraps_origin = wraps_origin,
    
    left_gene = left_gene,
    
    right_gene = right_gene
    
  )
  
  if(!is_single_string(interval$left_gene)){
    stop("left_gene must be a single character string.")
  }
  
  if(!is_single_string(interval$right_gene)){
    stop("right_gene must be a single character string.")
  }
  
  class(interval) <- "Interval"
  
  return(
    validate_interval(interval)
  )
  
}


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Assert Interval class
#
# Description:
#   Checks that an object inherits from the Interval class.
#
# Arguments:
#   interval : Object to test.
#
# Returns:
#   Invisibly returns the Interval object.
#--------------------------------------------------------------

assert_interval <- function(interval){
  
  if(!inherits(interval, "Interval")){
    
    stop(
      "Input must be an Interval object."
    )
    
  }
  
  invisible(interval)
  
}


#--------------------------------------------------------------
# Validate an Interval object
#
# Description:
#   Performs integrity checks on an Interval object.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   The validated Interval object.
#--------------------------------------------------------------

validate_interval <- function(interval){
  
  assert_interval(interval)
  
  required <- c(
    
    "start",
    
    "end",
    
    "width",
    
    "wraps_origin",
    
    "left_gene",
    
    "right_gene"
    
  )
  
  missing <- required[
    !required %in% names(interval)
  ]
  
  if(length(missing) > 0){
    
    stop(
      
      "Missing interval fields: ",
      
      paste(
        missing,
        collapse = ", "
      )
      
    )
    
  }
  
  if(
    !is_valid_coordinate(interval$start) ||
    !is_valid_coordinate(interval$end)
  ){
    stop("Interval coordinates are invalid.")
  }
  
  if(!is.logical(interval$wraps_origin) ||
     length(interval$wraps_origin) != 1){
    
    stop(
      "wraps_origin must be TRUE or FALSE."
    )
    
  }
  
  if(interval$start < 1){
    
    stop(
      "Interval starts before coordinate 1."
    )
    
  }
  
  if(!interval$wraps_origin &&
     interval$end < interval$start){
    
    stop(
      "Interval end precedes interval start."
    )
    
  }
  
  if(!interval$wraps_origin){
    
    expected_width <-
      interval$end -
      interval$start +
      1
    
    if(interval$width != expected_width){
      
      stop(
        "Interval width is inconsistent with coordinates."
      )
      
    }
    
  }
  
  return(interval)
  
}


#==============================================================
# Accessors
#==============================================================

#--------------------------------------------------------------
# Retrieve interval start coordinate
#
# Description:
#   Returns the first coordinate of an Interval.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   Integer start coordinate.
#--------------------------------------------------------------

get_interval_start <- function(interval){
  
  assert_interval(interval)
  
  interval$start
  
}


#--------------------------------------------------------------
# Retrieve interval end coordinate
#
# Description:
#   Returns the last coordinate of an Interval.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   Integer end coordinate.
#--------------------------------------------------------------

get_interval_end <- function(interval){
  
  assert_interval(interval)
  
  interval$end
  
}


#--------------------------------------------------------------
# Retrieve interval width
#
# Description:
#   Returns the width of an Interval.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   Interval width in base pairs.
#--------------------------------------------------------------

get_interval_width <- function(interval){
  
  assert_interval(interval)
  
  interval$width
  
}


#--------------------------------------------------------------
# Retrieve the left flanking gene
#
# Description:
#   Returns the name of the left flanking gene.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_left_gene <- function(interval){
  
  assert_interval(interval)
  
  interval$left_gene
  
}


#--------------------------------------------------------------
# Retrieve the right flanking gene
#
# Description:
#   Returns the name of the right flanking gene.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   Character string.
#--------------------------------------------------------------

get_right_gene <- function(interval){
  
  assert_interval(interval)
  
  interval$right_gene
  
}


#--------------------------------------------------------------
# Determine whether an Interval wraps the genome origin
#
# Description:
#   Indicates whether an Interval crosses the origin of a
#   circular genome.
#
# Arguments:
#   interval : An Interval object.
#
# Returns:
#   TRUE or FALSE.
#--------------------------------------------------------------

interval_wraps_origin <- function(interval){
  
  assert_interval(interval)
  
  interval$wraps_origin
  
}


#==============================================================
# Methods
#==============================================================

#--------------------------------------------------------------
# Print an Interval object
#
# Description:
#   Displays a concise summary of an Interval object.
#
# Arguments:
#   x : An Interval object.
#
# Returns:
#   Invisibly returns the Interval object.
#--------------------------------------------------------------

print.Interval <- function(x, ...){
  
  x <- validate_interval(x)
  
  cat("\n")
  
  cat("Interval Object\n")
  
  cat("---------------------------\n")
  
  cat("Start        :", x$start, "\n")
  
  cat("End          :", x$end, "\n")
  
  cat("Width        :", x$width, "bp\n")
  
  cat("Wraps origin :", x$wraps_origin, "\n")
  
  cat("Left gene    :", x$left_gene, "\n")
  
  cat("Right gene   :", x$right_gene, "\n")
  
  invisible(x)
  
}