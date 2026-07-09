###############################################################
##
## report.R
##
## Version : 1.0
## Status  : Development
##
## Purpose:
##   Generate human-readable reports for barcode analyses.
##
## Responsibilities:
##   - Construct report objects
##   - Generate QC reports
##   - Generate analysis summary reports
##   - Export reports as TXT, Markdown or HTML
##
###############################################################


#==============================================================
# Constructors
#==============================================================

#--------------------------------------------------------------
# Create a report object
#
# Description:
#   Constructs a cpBarcodeReport object containing formatted
#   report text and metadata.
#
# Arguments:
#   title   : Report title.
#   content : Character vector containing report text.
#
# Returns:
#   A validated cpBarcodeReport object.
#--------------------------------------------------------------

new_report <- function(title,
                       content){
  
  if(!is_single_string(title)){
    
    stop(
      "'title' must be a single character string."
    )
    
  }
  
  if(!is.character(content)){
    
    stop(
      "'content' must be a character vector."
    )
    
  }
  
  report <- list(
    
    title = title,
    
    content = content,
    
    created = Sys.time(),
    
    package = PACKAGE_NAME
    
  )
  
  class(report) <- "cpBarcodeReport"
  
  return(report)
  
}


#==============================================================
# QC Reports
#==============================================================

#--------------------------------------------------------------
# Generate a barcode quality-control report
#
# Description:
#   Converts a cpBarcodeQC object into a formatted report.
#
# Arguments:
#   qc : A cpBarcodeQC object.
#
# Returns:
#   A cpBarcodeReport object.
#--------------------------------------------------------------

report_qc <- function(qc){
  
  assert_cpBarcodeQC(qc)
  
  warnings <- if(length(qc$warnings) == 0){
    
    "None"
    
  } else{
    
    paste(
      qc$warnings,
      collapse = "\n"
    )
    
  }
  
  content <- c(
    
    paste("Marker              :", qc$marker),
    
    paste("Organism            :", qc$organism),
    
    paste("Accession           :", qc$accession),
    
    "",
    
    paste("Sequence length     :", qc$length,
          "bp"),
    
    sprintf(
      "GC content          : %.2f %%",
      qc$gc_content
    ),
    
    sprintf(
      "AT content          : %.2f %%",
      qc$at_content
    ),
    
    paste(
      "Ambiguous bases     :",
      qc$ambiguous
    ),
    
    sprintf(
      "Ambiguous (%%)       : %.2f %%",
      qc$ambiguous_percentage
    ),
    
    paste(
      "Longest homopolymer :",
      qc$longest_homopolymer,
      "bp"
    ),
    
    "",
    
    "Warnings",
    
    "--------",
    
    warnings
    
  )
  
  new_report(
    
    title = "DNA Barcode QC Report",
    
    content = content
    
  )
  
}


#==============================================================
# Summary Reports
#==============================================================

#--------------------------------------------------------------
# Generate an analysis summary report
#
# Description:
#   Converts a summary object into a formatted report.
#
# Arguments:
#   summary : Summary object.
#
# Returns:
#   A cpBarcodeReport object.
#--------------------------------------------------------------

report_summary <- function(summary){
  
  content <- capture.output(
    
    print(summary)
    
  )
  
  new_report(
    
    title = "Barcode Analysis Summary",
    
    content = content
    
  )
  
}


#==============================================================
# Export
#==============================================================

#--------------------------------------------------------------
# Write a report to disk
#
# Description:
#   Exports a cpBarcodeReport object as a plain text,
#   Markdown or HTML document.
#
# Arguments:
#   report : A cpBarcodeReport object.
#   file   : Output filename.
#   format : Optional output format. If NULL, the format is
#            inferred from the filename extension.
#
# Returns:
#   Invisibly returns the output filename.
#--------------------------------------------------------------

write_report <- function(report,
                         file,
                         format = NULL){
  
  if(!inherits(report,
               "cpBarcodeReport")){
    
    stop(
      "Input must be a cpBarcodeReport object."
    )
    
  }
  
  if(!is_single_string(file)){
    stop(
      "'file' must be a single character string."
    )
  }
  
  if(is.null(format)){
    
    format <- tolower(
      
      tools::file_ext(file)
      
    )
    
  }
  
  format <- sub("^\\.",
                "",
                format)
  
  if(!format %in%
     c("txt",
       "md",
       "html")){
    
    stop(
      "Unsupported report format."
    )
    
  }
  
  ensure_directory(
    
    dirname(file)
    
  )
  
  output <- switch(
    
    format,
    
    txt = c(
      
      report$title,
      
      paste(
        
        rep("=",
            nchar(report$title)),
        
        collapse = ""
        
      ),
      
      "",
      
      report$content
      
    ),
    
    md = c(
      
      paste0("# ",
             report$title),
      
      "",
      
      paste0(
        "**Generated:** ",
        format(
          report$created,
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      
      "",
      
      report$content
      
    ),
    
    html = c(
      
      "<!DOCTYPE html>",
      
      "<html>",
      
      "<head>",
      
      "<meta charset=\"UTF-8\">",
      
      paste0(
        "<title>",
        report$title,
        "</title>"
      ),
      
      "</head>",
      
      "<body>",
      
      paste0(
        "<h1>",
        report$title,
        "</h1>"
      ),
      
      paste0(
        
        "<p><strong>Generated:</strong> ",
        
        format(
          report$created,
          "%Y-%m-%d %H:%M:%S"
        ),
        
        "</p>"
        
      ),
      
      "<pre>",
      
      paste(
        
        report$content,
        
        collapse = "\n"
        
      ),
      
      "</pre>",
      
      "</body>",
      
      "</html>"
      
    )
    
  )
  
  writeLines(
    
    output,
    
    con = file
    
  )
  
  invisible(file)
  
}


#==============================================================
# Methods
#==============================================================

#--------------------------------------------------------------
# Print a report
#
# Description:
#   Displays a formatted report.
#
# Arguments:
#   x : A cpBarcodeReport object.
#   ... : Additional arguments (unused).
#
# Returns:
#   Invisibly returns the report.
#--------------------------------------------------------------

print.cpBarcodeReport <- function(x,
                                  ...){
  
  assert_cpBarcodeReport <- function(report){

  if(!inherits(report, "cpBarcodeReport")){

    stop(
      "Input must be a cpBarcodeReport object."
    )

  }

  invisible(report)

}
  
  assert_cpBarcodeReport(x)
  
  cat("\n")
  
  cat(x$title,
      "\n")
  
  cat(
    
    paste(
      
      rep("=",
          nchar(x$title)),
      
      collapse = ""
      
    ),
    
    "\n\n",
    
    sep = ""
    
  )
  
  cat(
    
    paste(
      
      x$content,
      
      collapse = "\n"
      
    ),
    
    "\n"
    
  )
  
  invisible(x)
  
}


