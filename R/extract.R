###############################################################
##
## extract.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Extract DNA barcode markers from chloroplast genome
##   objects.
##
## Responsibilities:
##   - Compute genomic intervals
##   - Extract DNA subsequences
##   - Extract gene markers
##   - Extract intergenic markers
##   - Dispatch marker extraction workflows
##
###############################################################

#==============================================================
# Interval Utilities
#==============================================================

#--------------------------------------------------------------
# Compute the interval between two genomic features
#
# Description:
#   Computes the genomic interval separating two annotated
#   features and returns it as an Interval object.
#
# Arguments:
#   feature1 : First parsed feature.
#   feature2 : Second parsed feature.
#
# Returns:
#   A validated Interval object.
#
# Notes:
#   Joined features (e.g. rps12) cannot currently be used for
#   interval calculations because they have undefined start/end
#   coordinates.
#--------------------------------------------------------------

compute_interval <- function(feature1,
                             feature2){
  
  #---------------------------------------------
  # Validate inputs
  #---------------------------------------------
  
  required <- c(
    "start",
    "end",
    "gene"
  )
  
  for(feature in list(feature1, feature2)){
    
    if(!is.list(feature)){
      
      stop(
        "Parsed feature must be a list."
      )
      
    }
    
    missing <- setdiff(
      required,
      names(feature)
    )
    
    if(length(missing) > 0){
      
      stop(
        "Parsed feature is missing required fields: ",
        paste(missing, collapse = ", ")
      )
      
    }
    
  }
  
  if(
    !is_valid_coordinate(feature1$start) ||
    !is_valid_coordinate(feature1$end) ||
    !is_valid_coordinate(feature2$start) ||
    !is_valid_coordinate(feature2$end)
  ){
    
    stop(
      "Cannot compute interval from joined or incomplete features."
    )
    
  }
  
  #---------------------------------------------
  # Determine which feature occurs first
  #---------------------------------------------
  
  if (feature1$start <= feature2$start) {
    
    left <- feature1
    
    right <- feature2
    
  } else {
    
    left <- feature2
    
    right <- feature1
    
  }
  
  interval_start <- left$end + 1
  
  interval_end <- right$start - 1
  
  wraps_origin <- interval_start > interval_end
  
  return(
    
    new_interval(
      
      start = interval_start,
      
      end = interval_end,
      
      wraps_origin = wraps_origin,
      
      left_gene = left$gene,
      
      right_gene = right$gene
      
    )
    
  )
  
}

#==============================================================
# Sequence Utilities
#==============================================================

#--------------------------------------------------------------
# Extract a DNA subsequence
#
# Description:
#   Extracts a DNA subsequence from a genome sequence using a
#   validated Interval object. Circular genomes that span the
#   origin are handled automatically.
#
# Arguments:
#   sequence : Genome sequence.
#   interval : An Interval object.
#
# Returns:
#   Character string containing the extracted DNA sequence.
#--------------------------------------------------------------

extract_subsequence <- function(sequence,
                                interval){
  
  assert_interval(interval)
  
  if(!is_valid_dna(sequence)){
    
    stop(
      "'sequence' must be a valid DNA sequence."
    )
    
  }
  
  if (!interval$wraps_origin) {
    
    dna <- substr(
      sequence,
      interval$start,
      interval$end
    )
    
    return(dna)
    
  }
  
  # Circular genome
  left_part <- substr(
    sequence,
    interval$start,
    nchar(sequence)
  )
  
  right_part <- substr(
    sequence,
    1,
    interval$end
  )
  
  return(
    paste0(
      left_part,
      right_part
    )
  )
  
}


#==============================================================
# Marker Extraction
#==============================================================

#--------------------------------------------------------------
# Extract a gene barcode marker from a cpGenome object
#
# Description:
#   Locates a gene annotation, extracts the corresponding DNA
#   sequence, and returns a cpBarcode object.
#
# Arguments:
#   genome      : A cpGenome object.
#   marker_info : Marker definition from MARKER_DICTIONARY.
#   options     : Pipeline options.
#
# Returns:
#   A validated cpBarcode object.
#--------------------------------------------------------------

extract_gene <- function(genome,
                         marker_info,
                         options){
  
  assert_cpGenome(genome)
  
  options <- validate_pipeline_options(options)
  
  #---------------------------------------------
  # Locate the gene
  #---------------------------------------------
  
  feature <- get_gene_feature(
    
    feature_table = get_feature_table(genome),
    
    gene_name = marker_info$gene
    
  )
  
  #---------------------------------------------
  # Extract the DNA sequence
  #---------------------------------------------
  
  dna <- substr(
    
    get_sequence(genome),
    
    feature$start,
    
    feature$end
    
  )
  
  #---------------------------------------------
  # Reverse complement if requested
  #---------------------------------------------
  
  if(options$reverse_complement &&
     feature$strand == "-"){
    
    dna <- reverse_complement(dna)
    
  }
  
  #---------------------------------------------
  # Construct barcode object
  #---------------------------------------------
  
  barcode <- new_cpBarcode(
    
    marker = marker_info$marker,
    
    accession = get_accession(genome),
    
    version = get_version(genome),
    
    organism = get_organism(genome),
    
    feature_type = feature$type,
    
    sequence = dna,
    
    coordinates = list(
      
      start = feature$start,
      
      end = feature$end
      
    ),
    
    strand = feature$strand,
    
    metadata = list(
      
      gene = marker_info$gene,
      
      product = feature$product,
      
      locus_tag = feature$locus_tag,
      
      protein_id = feature$protein_id
      
    )
    
  )
  
  validate_barcode(barcode)
  
}


#--------------------------------------------------------------
# Extract an intergenic barcode marker from a cpGenome object
#
# Description:
#   Locates two flanking genes, computes the intergenic interval,
#   extracts the DNA sequence, and returns a cpBarcode object.
#
# Arguments:
#   genome      : cpGenome object.
#   marker_info : Marker definition from MARKER_DICTIONARY.
#   options     : Pipeline options.
#
# Returns:
#   A cpBarcode object.
#--------------------------------------------------------------

extract_intergenic <- function(genome,
                               marker_info,
                               options){
  
  assert_cpGenome(genome)
  
  options <- validate_pipeline_options(options)
  
  #---------------------------------------------
  # Marker information
  #---------------------------------------------
  
  feature_type <- FEATURE_INTERGENIC
  
  #---------------------------------------------
  # Locate the flanking genes
  #---------------------------------------------
  
  left_feature <- get_gene_feature(
    
    feature_table = get_feature_table(genome),
    
    gene_name = marker_info$left_gene
    
  )
  
  right_feature <- get_gene_feature(
    
    feature_table = get_feature_table(genome),
    
    gene_name = marker_info$right_gene
    
  )
  
  #---------------------------------------------
  # Compute the intergenic interval
  #---------------------------------------------
  
  interval <- compute_interval(
    
    feature1 = left_feature,
    
    feature2 = right_feature
    
  )
  
  #---------------------------------------------
  # Extract the DNA sequence
  #---------------------------------------------
  
  dna <- extract_subsequence(
    
    sequence = get_sequence(genome),
    
    interval = interval
    
  )
  
  #---------------------------------------------
  # Construct the barcode object
  #---------------------------------------------
  
  barcode <- new_cpBarcode(
    
    marker = marker_info$marker,
    
    accession = get_accession(genome),
    
    version = get_version(genome),
    
    organism = get_organism(genome),
    
    feature_type = feature_type,
    
    sequence = dna,
    
    coordinates = list(
      
      start = get_interval_start(interval),
      
      end = get_interval_end(interval)
      
    ),
    
    strand = STRAND_FORWARD,
    
    metadata = list(
      
      left_gene = marker_info$left_gene,
      
      right_gene = marker_info$right_gene,
      
      wraps_origin = interval_wraps_origin(interval)
      
    )
    
  )
  
  return(
    
    validate_barcode(barcode)
    
  )
  
}


#==============================================================
# High-level API
#==============================================================

#--------------------------------------------------------------
# Extract a DNA barcode marker
#
# Description:
#   Dispatches extraction to the appropriate method according
#   to the marker definition.
#
# Arguments:
#   genome  : A cpGenome object.
#   marker  : Marker name.
#   options : Pipeline options.
#
# Returns:
#   Marker extraction result.
#--------------------------------------------------------------

extract_marker <- function(genome,
                           marker,
                           options = PIPELINE_OPTIONS){
  
  #---------------------------------------------
  # Validate input
  #---------------------------------------------
  
  assert_cpGenome(genome)
  
  options <- validate_pipeline_options(options)
  
  #---------------------------------------------
  # Retrieve marker definition
  #---------------------------------------------
  
  marker_info <- MARKER_DICTIONARY[[marker]]
  
  if (is.null(marker_info)) {
    
    stop(
      "Unknown marker: ",
      marker
    )
    
  }
  
  if (is.null(marker_info$class)) {
    
    stop(
      "Marker definition has no class."
    )
    
  }
  
  #---------------------------------------------
  # Dispatch to extraction method
  #---------------------------------------------
  
  if(marker_info$class == MARKER_CLASS_GENE){
    
    return(
      extract_gene(
        genome,
        marker_info,
        options
      )
    )
    
  }
  
  if(marker_info$class == MARKER_CLASS_INTERGENIC){
    
    return(
      extract_intergenic(
        genome,
        marker_info,
        options
      )
    )
    
  }
  
  if(marker_info$class == MARKER_CLASS_NUCLEAR){
    
    stop(
      "NUCLEAR markers are not implemented yet."
    )
    
  }
  
  stop(
    "Unsupported marker class."
  )

}