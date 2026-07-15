###############################################################
##
## batch.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Perform batch extraction of chloroplast DNA barcode
##   markers from multiple GenBank files.
##
## Responsibilities:
##   - Extract one marker from multiple genomes
##   - Extract multiple markers from multiple genomes
##   - Continue processing when individual files fail
##   - Optionally export extracted barcodes
##
###############################################################

#==============================================================
# Batch Processing
#==============================================================

#--------------------------------------------------------------
# Extract one barcode marker from multiple genomes
#
# Description:
#   Iterates over a collection of GenBank files, extracting
#   the same barcode marker from each genome. Extraction
#   failures are reported as warnings without terminating the
#   batch process.
#
# Arguments:
#   files      : Character vector of GenBank filenames.
#   marker     : Barcode marker name.
#   output_dir : Optional output directory. If NULL, extracted
#                barcodes are returned but not exported.
#   options    : Pipeline options.
#
# Returns:
#   Named list of cpBarcode objects. Files that fail return
#   NULL.
#--------------------------------------------------------------

batch_extract_barcode <- function(files,
                                  marker,
                                  output_dir = NULL,
                                  format = "fasta",
                                  options = PIPELINE_OPTIONS){
  
  #---------------------------------------------
  # Validate inputs
  #---------------------------------------------
  
  if(!is.character(files) ||
     length(files) == 0 ||
     any(!vapply(files,
                 is_single_string,
                 logical(1)))){
    
    stop(
      "'files' must be a non-empty character vector of file paths."
    )
    
  }
  
  if(!is_single_string(marker)){
    
    stop(
      "'marker' must be a single character string."
    )
    
  }
  
  if(!is_single_string(format)){
    
    stop(
      "'format' must be a single character string."
    )
    
  }
  
  if(!is.null(output_dir)){
    
    if(!is_single_string(output_dir)){
      
      stop(
        "'output_dir' must be a single character string."
      )
      
    }
    
    ensure_directory(output_dir)
    
  }
  
  options <- validate_pipeline_options(options)
  
  results <- vector(
    mode = "list",
    length = length(files)
  )
  
  names(results) <- basename(files)
  
  for(i in seq_along(files)){
    
    file <- files[i]
    
    message(
      "[",
      i,
      "/",
      length(files),
      "] ",
      basename(file)
    )
    
    result <- tryCatch(
      
      {
        
        barcode <- extract_barcode(
          
          file = file,
          
          marker = marker,
          
          options = options
          
        )
        
        if(!is.null(output_dir)){
          
          outfile <- file.path(
            
            output_dir,
            
            paste0(
              tools::file_path_sans_ext(
                basename(file)
              ),
              ".fasta"
            )
            
          )
          
          write_barcode(
            
            barcode = barcode,
            
            file = outfile,
            
            format = format
            
          )
          
        }
        
        barcode
        
      },
      
      error = function(e){
        
        warning(
          
          basename(file),
          
          ": ",
          
          conditionMessage(e),
          
          call. = FALSE
          
        )
        
        NULL
        
      }
      
    )
    
    results[i] <- list(result)
    
  }
  
  return(results)
  
}

#--------------------------------------------------------------
# Extract multiple barcode markers from multiple genomes
#
# Description:
#   Iterates over multiple GenBank files, extracting every
#   requested barcode marker from each genome. Individual
#   extraction failures do not interrupt the remaining batch.
#
# Arguments:
#   files      : Character vector of GenBank filenames.
#   markers    : Character vector of marker names.
#   output_dir : Optional output directory.
#   options    : Pipeline options.
#
# Returns:
#   Named list where each element contains the extracted
#   barcodes for one genome.
#--------------------------------------------------------------

batch_extract_barcodes <- function(files,
                                   markers,
                                   output_dir = NULL,
                                   format = "fasta",
                                   options = PIPELINE_OPTIONS){

  #---------------------------------------------
  # Validate inputs
  #---------------------------------------------

  if(!is.character(files) ||
     length(files) == 0 ||
     any(!vapply(files,
                 is_single_string,
                 logical(1)))){

    stop(
      "'files' must be a non-empty character vector of file paths."
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

  if(!is_single_string(format)){

    stop(
      "'format' must be a single character string."
    )

  }

  if(!is.null(output_dir)){

    if(!is_single_string(output_dir)){

      stop(
        "'output_dir' must be a single character string."
      )

    }

    ensure_directory(output_dir)

  }

  options <- validate_pipeline_options(options)

  results <- vector(
    mode = "list",
    length = length(files)
  )

  names(results) <- basename(files)

  for(i in seq_along(files)){

    file <- files[i]

    message(

      "[",

      i,

      "/",

      length(files),

      "] ",

      basename(file)

    )

    result <- tryCatch(

      {

        barcodes <- extract_barcodes(

          file = file,

          markers = markers,

          output_dir = NULL,

          format = format,

          options = options

        )

        barcodes

      },

      error = function(e){

        warning(

          basename(file),

          ": ",

          conditionMessage(e),

          call. = FALSE

        )

        NULL

      }

    )

    results[i] <- list(result)

  }

  #---------------------------------------------
  # Export combined marker files
  #---------------------------------------------

  if(!is.null(output_dir)){

    marker_names <- unique(

      unlist(

        lapply(

          results,

          names

        ),

        use.names = FALSE

      )

    )

    for(marker_name in marker_names){

      marker_barcodes <- lapply(

        results,

        function(result){

          if(is.null(result)){

            return(NULL)

          }

          result[[marker_name]]

        }

      )

      marker_barcodes <- Filter(

        Negate(is.null),

        marker_barcodes

      )

      if(length(marker_barcodes) == 0){

        next

      }

      output_file <- file.path(

        output_dir,

        paste0(

          marker_name,

          ".",

          tolower(format)

        )

      )

      if(tolower(format) == "fasta"){

        fasta_lines <- unlist(

          lapply(

            marker_barcodes,

            as_fasta

          ),

          use.names = FALSE

        )

        writeLines(

          fasta_lines,

          output_file

        )

      }else if(tolower(format) == "csv"){

        barcode_data <- do.call(

          rbind,

          lapply(

            marker_barcodes,

            as_dataframe

          )

        )

        utils::write.csv(

          barcode_data,

          output_file,

          row.names = FALSE

        )

      }else{

        stop(

          "Combined batch export currently supports only 'fasta' and 'csv'."

        )

      }

    }

  }

  return(results)

}
