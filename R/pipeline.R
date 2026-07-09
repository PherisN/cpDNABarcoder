###############################################################
##
## pipeline.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Provide high-level user interfaces for DNA barcode
##   extraction workflows.
##
## Responsibilities:
##   - Read chloroplast genomes
##   - Extract a single DNA barcode
##   - Extract multiple DNA barcodes
##   - Export extracted barcode sequences
##
###############################################################

#==============================================================
# High-level Extraction
#==============================================================

#--------------------------------------------------------------
# Extract a single DNA barcode
#
# Description:
#   Reads a GenBank file, extracts the requested DNA barcode
#   marker, optionally exports it, and returns a validated
#   cpBarcode object.
#
# Arguments:
#   file    : Path to a GenBank file.
#   marker  : Marker name.
#   output  : Optional output filename.
#   format  : Output format.
#   options : Pipeline options.
#
# Returns:
#   A validated cpBarcode object.
#--------------------------------------------------------------

extract_barcode <- function(file,
                            marker,
                            output = NULL,
                            format = "fasta",
                            options = PIPELINE_OPTIONS){
  
  #---------------------------------------------
  # Validate inputs
  #---------------------------------------------
  
  options <- validate_pipeline_options(options)
  
  if(!is_single_string(marker)){
    
    stop(
      "'marker' must be a single character string."
    )
    
  }
  
  #---------------------------------------------
  # Read genome
  #---------------------------------------------
  
  genome <- read_cpGenome(file)
  
  #---------------------------------------------
  # Extract barcode
  #---------------------------------------------
  
  barcode <- validate_barcode(
    
    extract_marker(
      
      genome = genome,
      
      marker = marker,
      
      options = options
      
    )
    
  )
  
  #---------------------------------------------
  # Export if requested
  #---------------------------------------------
  
  if(!is.null(output)){
    
    write_barcode(
      
      barcode = barcode,
      
      file = output,
      
      format = format
      
    )
    
  }
  
  return(barcode)
  
}



#--------------------------------------------------------------
# Extract multiple DNA barcode markers
#
# Description:
#   Reads a GenBank file, extracts multiple DNA barcode
#   markers, optionally exports each barcode, and returns a
#   named list of validated cpBarcode objects.
#
# Arguments:
#   file       : Path to a GenBank file.
#   markers    : Character vector of marker names.
#   output_dir : Optional output directory.
#   format     : Output format.
#   options    : Pipeline options.
#
# Returns:
#   A named list of validated cpBarcode objects.
#--------------------------------------------------------------

extract_barcodes <- function(file,
                             markers,
                             output_dir = NULL,
                             format = "fasta",
                             options = PIPELINE_OPTIONS){
  
  #---------------------------------------------
  # Validate inputs
  #---------------------------------------------
  
  options <- validate_pipeline_options(options)
  
  if(!is.null(output_dir) &&
     !is_single_string(output_dir)){
    
    stop(
      "'output_dir' must be a single character string."
    )
    
  }
  
  if(!is.character(markers) ||
     length(markers) == 0 ||
     any(!vapply(markers,
                 is_single_string,
                 logical(1)))){
    
    stop(
      "'markers' must be a non-empty character vector of non-missing strings."
    )
    
  }
  
  #---------------------------------------------
  # Read genome
  #---------------------------------------------
  
  genome <- read_cpGenome(file)
  
  #---------------------------------------------
  # Prepare output
  #---------------------------------------------
  
  barcodes <- vector(
    
    mode = "list",
    
    length = length(markers)
    
  )
  
  names(barcodes) <- markers
  
  if(!is.null(output_dir)){
    
    ensure_directory(output_dir)
    
  }
  
  #---------------------------------------------
  # Extract markers
  #---------------------------------------------
  
  for(i in seq_along(markers)){
    
    marker <- markers[i]
    
    barcode <- validate_barcode(
      
      extract_marker(
        
        genome = genome,
        
        marker = marker,
        
        options = options
        
      )
      
    )
    
    if(!is.null(output_dir)){
      
      outfile <- file.path(
        
        output_dir,
        
        paste0(
          marker,
          ".",
          tolower(format)
        )
        
      )
      
      write_barcode(
        
        barcode = barcode,
        
        file = outfile,
        
        format = format
        
      )
      
    }
    
    barcodes[[i]] <- barcode
    
  }
  
  return(barcodes)
  
}