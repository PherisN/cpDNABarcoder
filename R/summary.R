###############################################################
##
## summary.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Summarize extracted DNA barcode objects and collections.
##
## Responsibilities:
##   - Summarize individual barcodes
##   - Summarize barcode collections
##   - Compute collection-level statistics
##   - Display concise summary reports
##
###############################################################

#==============================================================
# Barcode Summaries
#==============================================================

#--------------------------------------------------------------
# Summarize a cpBarcode object
#
# Description:
#   Computes a concise summary of a cpBarcode object by
#   combining barcode metadata with sequence quality metrics.
#
# Arguments:
#   object : A cpBarcode object.
#   ...    : Reserved for future extensions.
#
# Returns:
#   A cpBarcodeSummary object.
#--------------------------------------------------------------

summary.cpBarcode <- function(object, ...){
  
  barcode <- validate_barcode(object)
  
  qc <- qc_barcode(barcode)
  
  summary <- list(
    
    marker = get_marker(barcode),
    
    accession = get_accession(barcode),
    
    organism = get_organism(barcode),
    
    feature_type = get_feature_type(barcode),
    
    coordinates = get_coordinates(barcode),
    
    strand = get_strand(barcode),
    
    length = qc$length,
    
    gc_content = qc$gc_content,
    
    at_content = qc$at_content,
    
    ambiguous = qc$ambiguous,
    
    ambiguous_percentage =
      qc$ambiguous_percentage,
    
    longest_homopolymer =
      qc$longest_homopolymer
    
  )
  
  class(summary) <- "cpBarcodeSummary"
  
  summary
  
}


#--------------------------------------------------------------
# Print a barcode summary
#
# Description:
#   Displays a concise summary of a cpBarcodeSummary object.
#
# Arguments:
#   x : A cpBarcodeSummary object.
#
# Returns:
#   Invisibly returns the object.
#--------------------------------------------------------------

print.cpBarcodeSummary <- function(x, ...){
  
  cat("\n")
  
  cat("Barcode Summary\n")
  
  cat("---------------------------\n")
  
  cat("Marker       :", x$marker, "\n")
  
  cat("Organism     :", x$organism, "\n")
  
  cat("Accession    :", x$accession, "\n")
  
  cat("Feature type :", x$feature_type, "\n")
  
  cat(
    "Coordinates  :",
    paste(
      x$coordinates$start,
      x$coordinates$end,
      sep = "-"
    ),
    "\n"
  )
  
  cat("Strand       :", x$strand, "\n")
  
  cat("Length       :", x$length, "bp\n")
  
  cat(sprintf(
    "GC content   : %.2f %%\n",
    x$gc_content
  ))
  
  cat(sprintf(
    "AT content   : %.2f %%\n",
    x$at_content
  ))
  
  cat(
    "Ambiguous    :",
    x$ambiguous,
    "\n"
  )
  
  cat(sprintf(
    "Ambiguous %% : %.2f %%\n",
    x$ambiguous_percentage
  ))
  
  cat(
    "Homopolymer  :",
    x$longest_homopolymer,
    "bp\n"
  )
  
  invisible(x)
  
}


#==============================================================
# Collection Summaries
#==============================================================

#--------------------------------------------------------------
# Summarize a collection of barcodes
#
# Description:
#   Computes collection-level statistics for a list of
#   cpBarcode objects.
#
# Arguments:
#   barcodes : List of cpBarcode objects.
#
# Returns:
#   A cpBarcodeCollectionSummary object.
#--------------------------------------------------------------

summarize_barcodes <- function(barcodes){
  
  if(!is.list(barcodes)){
    
    stop(
      "barcodes must be supplied as a list."
    )
    
  }
  
  if(length(barcodes) == 0){
    
    stop(
      "'barcodes' cannot be empty."
    )
    
  }
  
  if(any(vapply(
    barcodes,
    function(x) !inherits(x, "cpBarcode"),
    logical(1)
  ))){
    
    stop(
      "All elements of 'barcodes' must be cpBarcode objects."
    )
    
  }
  
  barcodes <- Filter(Negate(is.null), barcodes)
  
  if(length(barcodes) == 0){
    
    stop(
      "No valid barcode objects supplied."
    )
    
  }
  
  summaries <- lapply(
    barcodes,
    summary
  )
  
  lengths <- vapply(
    summaries,
    function(x) x$length,
    numeric(1)
  )
  
  gc <- vapply(
    summaries,
    function(x) x$gc_content,
    numeric(1)
  )
  
  markers <- vapply(
    summaries,
    function(x) x$marker,
    character(1)
  )
  
  organisms <- vapply(
    summaries,
    function(x) x$organism,
    character(1)
  )
  
  ambiguous <- vapply(
    summaries,
    function(x) x$ambiguous,
    numeric(1)
  )
  
  report <- list(
    
    total_barcodes = length(barcodes),
    
    marker_frequency = table(markers),
    
    organism_count = length(unique(organisms)),
    
    mean_length = mean(lengths),
    
    min_length = min(lengths),
    
    max_length = max(lengths),
    
    mean_gc = mean(gc),
    
    total_ambiguous = sum(ambiguous)
    
  )
  
  class(report) <- "cpBarcodeCollectionSummary"
  
  report
  
}


#--------------------------------------------------------------
# Print a collection summary
#
# Description:
#   Displays summary statistics for a barcode collection.
#
# Arguments:
#   x : A cpBarcodeCollectionSummary object.
#
# Returns:
#   Invisibly returns the object.
#--------------------------------------------------------------

print.cpBarcodeCollectionSummary <- function(x, ...){
  
  cat("\n")
  
  cat("Barcode Collection Summary\n")
  
  cat("---------------------------\n")
  
  cat(
    "Total barcodes :",
    x$total_barcodes,
    "\n"
  )
  
  cat(
    "Organisms      :",
    x$organism_count,
    "\n"
  )
  
  cat(sprintf(
    "Mean length    : %.1f bp\n",
    x$mean_length
  ))
  
  cat(
    "Minimum length :",
    x$min_length,
    "bp\n"
  )
  
  cat(
    "Maximum length :",
    x$max_length,
    "bp\n"
  )
  
  cat(sprintf(
    "Mean GC        : %.2f %%\n",
    x$mean_gc
  ))
  
  cat(
    "Ambiguous bases:",
    x$total_ambiguous,
    "\n"
  )
  
  cat("\nMarker frequencies\n")
  
  cat("---------------------------\n")
  
  print(x$marker_frequency)
  
  invisible(x)
  
}