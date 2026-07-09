###############################################################
##
## qc.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Compute sequence quality metrics for DNA barcode objects.
##
## Responsibilities:
##   - Compute GC and AT content
##   - Count ambiguous nucleotides
##   - Detect homopolymers
##   - Perform barcode quality assessment
##
###############################################################

#==============================================================
# Sequence Quality Metrics
#==============================================================

#--------------------------------------------------------------
# Check whether a DNA sequence contains only valid symbols
#
# Description:
#   Verifies that a DNA sequence contains only characters from
#   the supported IUPAC nucleotide alphabet.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   TRUE if valid, FALSE otherwise.
#--------------------------------------------------------------

is_valid_dna <- function(sequence){
  
  if(!is.character(sequence) ||
     length(sequence) != 1){
    
    return(FALSE)
    
  }
  
  sequence <- toupper(sequence)
  
  return(
    
    grepl(
      paste0(
        "^[",
        paste(IUPAC_BASES, collapse = ""),
        "]+$"
      ),
      sequence
    )
    
  )
  
}


#--------------------------------------------------------------
# Count ambiguous nucleotide symbols
#
# Description:
#   Counts nucleotide symbols that are not A, C, G or T.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Integer.
#--------------------------------------------------------------

count_ambiguous_bases <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    stop(
      "Sequence contains invalid DNA symbols."
    )
  }
  
  sequence <- toupper(sequence)
  
  bases <- strsplit(sequence, "")[[1]]
  
  sum(
    !bases %in% DNA_BASES
  )
  
}


#--------------------------------------------------------------
# Calculate nucleotide composition
#
# Description:
#   Counts each nucleotide symbol occurring in a DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Named integer vector.
#--------------------------------------------------------------

nucleotide_composition <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  bases <- strsplit(
    toupper(sequence),
    ""
  )[[1]]
  
  table(
    factor(
      bases,
      levels = IUPAC_BASES
    )
  )
  
}


#--------------------------------------------------------------
# Calculate GC content
#
# Description:
#   Calculates the GC percentage of a DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Numeric percentage.
#--------------------------------------------------------------

gc_content <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  bases <- strsplit(
    toupper(sequence),
    ""
  )[[1]]
  
  gc <- sum(
    bases %in% c("G", "C")
  )
  
  100 * gc / length(bases)
  
}


#--------------------------------------------------------------
# Calculate AT content
#
# Description:
#   Calculates the AT percentage of a DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Numeric percentage.
#--------------------------------------------------------------

at_content <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  bases <- strsplit(
    toupper(sequence),
    ""
  )[[1]]
  
  at <- sum(
    bases %in% c("A", "T")
  )
  
  100 * at / length(bases)
  
}


#--------------------------------------------------------------
# Calculate ambiguous base percentage
#
# Description:
#   Calculates the proportion of ambiguous nucleotide symbols.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Numeric percentage.
#--------------------------------------------------------------

ambiguous_percentage <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  100 *
    count_ambiguous_bases(sequence) /
    dna_length(sequence)
  
}


#--------------------------------------------------------------
# Determine the longest homopolymer run
#
# Description:
#   Determines the length of the longest consecutive run of
#   identical nucleotides.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Integer.
#--------------------------------------------------------------

longest_homopolymer <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  runs <- rle(
    strsplit(
      toupper(sequence),
      ""
    )[[1]]
  )
  
  max(runs$lengths)
  
}


#--------------------------------------------------------------
# Return sequence length
#
# Description:
#   Returns the number of nucleotides in a DNA sequence.
#
# Arguments:
#   sequence : DNA sequence.
#
# Returns:
#   Integer.
#--------------------------------------------------------------

sequence_length <- function(sequence){
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "Sequence contains invalid DNA symbols."
    )
    
  }
  
  dna_length(sequence)
  
}

#==============================================================
# Barcode QC
#==============================================================

#--------------------------------------------------------------
# Perform quality assessment of a DNA barcode
#
# Description:
#   Computes sequence quality metrics for a cpBarcode object.
#   The resulting quality report is intended for downstream
#   inspection and summary reporting.
#
# Arguments:
#   barcode : A cpBarcode object.
#
# Returns:
#   A named list containing:
#     - length
#     - gc_content
#     - at_content
#     - ambiguous
#     - ambiguous_percentage
#     - longest_homopolymer
#--------------------------------------------------------------

qc_barcode <- function(barcode){
  
  barcode <- validate_barcode(barcode)
  
  sequence <- get_sequence(barcode)
  
  report <- list(
    
    marker = get_marker(barcode),
    
    accession = get_accession(barcode),
    
    organism = get_organism(barcode),
    
    valid = is_valid_dna(sequence),
    
    length = sequence_length(sequence),
    
    gc_content = gc_content(sequence),
    
    at_content = at_content(sequence),
    
    ambiguous = count_ambiguous_bases(sequence),
    
    ambiguous_percentage =
      ambiguous_percentage(sequence),
    
    longest_homopolymer =
      longest_homopolymer(sequence),
    
    composition =
      nucleotide_composition(sequence),
    
    warnings = character(0)
    
  )
  
  if(report$ambiguous_percentage > 1){
    
    report$warnings <- c(
      
      report$warnings,
      
      "High proportion of ambiguous bases."
      
    )
    
  }
  
  if(report$gc_content < 20 ||
     report$gc_content > 50){
    
    report$warnings <- c(
      
      report$warnings,
      
      "GC content is outside the expected plastid range."
      
    )
    
  }
  
  class(report) <- "cpBarcodeQC"
  
  return(report)
  
}


#==============================================================
# Validation
#==============================================================

#--------------------------------------------------------------
# Validate a cpBarcodeQC object
#
# Description:
#   Ensures that an object inherits from the cpBarcodeQC class.
#
# Arguments:
#   report : Object to validate.
#
# Returns:
#   Invisibly returns the validated object.
#--------------------------------------------------------------

assert_cpBarcodeQC <- function(report){
  
  if(!inherits(report, "cpBarcodeQC")){
    
    stop(
      "Input must be a cpBarcodeQC object."
    )
    
  }
  
  invisible(report)
  
}


#==============================================================
# Methods
#==============================================================

#--------------------------------------------------------------
# Print a barcode QC report
#
# Description:
#   Displays a concise summary of quality-control metrics.
#
# Arguments:
#   x : A cpBarcodeQC object.
#
# Returns:
#   Invisibly returns the object.
#--------------------------------------------------------------

print.cpBarcodeQC <- function(x, ...){
  
  assert_cpBarcodeQC(x)
  
  cat("\n")
  
  cat("Barcode QC Report\n")
  
  cat("---------------------------\n")
  
  cat("Marker              :", x$marker, "\n")
  
  cat("Organism            :", x$organism, "\n")
  
  cat("Accession           :", x$accession, "\n")
  
  cat("Valid DNA           :", x$valid, "\n")
  
  cat("Length              :", x$length, "bp\n")
  
  cat(sprintf(
    "GC content          : %.2f %%\n",
    x$gc_content
  ))
  
  cat(sprintf(
    "AT content          : %.2f %%\n",
    x$at_content
  ))
  
  cat(
    "Ambiguous bases     :",
    x$ambiguous,
    "\n"
  )
  
  cat(sprintf(
    "Ambiguous (%%)       : %.2f %%\n",
    x$ambiguous_percentage
  ))
  
  cat(
    "Longest homopolymer :",
    x$longest_homopolymer,
    "bp\n"
  )
  
  cat("\nComposition\n")
  
  cat("---------------------------\n")
  
  print(x$composition)
  
  cat("\nWarnings\n")
  
  cat("---------------------------\n")
  
  if(length(x$warnings) == 0){
    
    cat("None\n")
    
  } else {
    
    cat(
      paste(
        x$warnings,
        collapse = "\n"
      ),
      "\n"
    )
    
  }
  
  invisible(x)
  
}