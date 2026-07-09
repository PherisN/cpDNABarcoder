make_test_report_barcode <- function(){
  
  valid_marker <- MARKER_DICTIONARY[[1]]$marker
  
  new_cpBarcode(
    marker = valid_marker,
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    feature_type = FEATURE_GENE,
    sequence = "AATTCCGGNN",
    coordinates = list(
      start = 1L,
      end = 10L
    ),
    strand = STRAND_FORWARD,
    metadata = list()
  )
  
}


testthat::test_that(
  "new_report constructs a cpBarcodeReport object",
  {
    
    report <- new_report(
      title = "Test Report",
      content = c(
        "Line 1",
        "Line 2"
      )
    )
    
    testthat::expect_true(
      inherits(
        report,
        "cpBarcodeReport"
      )
    )
    
    testthat::expect_equal(
      report$title,
      "Test Report"
    )
    
    testthat::expect_equal(
      report$content,
      c(
        "Line 1",
        "Line 2"
      )
    )
    
    testthat::expect_true(
      inherits(
        report$created,
        "POSIXct"
      )
    )
    
    testthat::expect_equal(
      report$package,
      PACKAGE_NAME
    )
    
  }
)


testthat::test_that(
  "new_report rejects invalid inputs",
  {
    
    testthat::expect_error(
      new_report(
        title = c("A", "B"),
        content = "Content"
      )
    )
    
    testthat::expect_error(
      new_report(
        title = NA_character_,
        content = "Content"
      )
    )
    
    testthat::expect_error(
      new_report(
        title = "Test Report",
        content = list("Content")
      )
    )
    
  }
)


testthat::test_that(
  "report_qc converts a cpBarcodeQC object to a report",
  {
    
    barcode <- make_test_report_barcode()
    
    qc <- qc_barcode(
      barcode
    )
    
    report <- report_qc(
      qc
    )
    
    testthat::expect_true(
      inherits(
        report,
        "cpBarcodeReport"
      )
    )
    
    testthat::expect_equal(
      report$title,
      "DNA Barcode QC Report"
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Marker",
          report$content
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Test organism",
          report$content
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "TEST001",
          report$content
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Warnings",
          report$content
        )
      )
    )
    
  }
)


testthat::test_that(
  "report_qc rejects non-QC objects",
  {
    
    testthat::expect_error(
      report_qc(
        list()
      )
    )
    
  }
)


testthat::test_that(
  "report_summary converts a barcode summary to a report",
  {
    
    barcode <- make_test_report_barcode()
    
    barcode_summary <- summary(
      barcode
    )
    
    report <- report_summary(
      barcode_summary
    )
    
    testthat::expect_true(
      inherits(
        report,
        "cpBarcodeReport"
      )
    )
    
    testthat::expect_equal(
      report$title,
      "Barcode Analysis Summary"
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Barcode Summary",
          report$content
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Test organism",
          report$content
        )
      )
    )
    
  }
)


testthat::test_that(
  "write_report writes a text report",
  {
    
    report <- new_report(
      title = "Test Report",
      content = c(
        "Line 1",
        "Line 2"
      )
    )
    
    file <- tempfile(
      fileext = ".txt"
    )
    
    result <- write_report(
      report,
      file
    )
    
    testthat::expect_equal(
      result,
      file
    )
    
    testthat::expect_true(
      file.exists(file)
    )
    
    lines <- readLines(
      file
    )
    
    testthat::expect_equal(
      lines[1],
      "Test Report"
    )
    
    testthat::expect_true(
      "Line 1" %in% lines
    )
    
    testthat::expect_true(
      "Line 2" %in% lines
    )
    
  }
)


testthat::test_that(
  "write_report writes a Markdown report",
  {
    
    report <- new_report(
      title = "Test Report",
      content = "Markdown content"
    )
    
    file <- tempfile(
      fileext = ".md"
    )
    
    write_report(
      report,
      file
    )
    
    lines <- readLines(
      file
    )
    
    testthat::expect_equal(
      lines[1],
      "# Test Report"
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Generated",
          lines
        )
      )
    )
    
    testthat::expect_true(
      "Markdown content" %in% lines
    )
    
  }
)


testthat::test_that(
  "write_report writes an HTML report",
  {
    
    report <- new_report(
      title = "Test Report",
      content = "HTML content"
    )
    
    file <- tempfile(
      fileext = ".html"
    )
    
    write_report(
      report,
      file
    )
    
    lines <- readLines(
      file
    )
    
    testthat::expect_true(
      "<!DOCTYPE html>" %in% lines
    )
    
    testthat::expect_true(
      "<h1>Test Report</h1>" %in% lines
    )
    
    testthat::expect_true(
      grepl(
        "HTML content",
        paste(
          lines,
          collapse = "\n"
        ),
        fixed = TRUE
      )
    )
    
  }
)


testthat::test_that(
  "write_report respects explicit format",
  {
    
    report <- new_report(
      title = "Test Report",
      content = "Explicit format content"
    )
    
    file <- tempfile(
      fileext = ".out"
    )
    
    write_report(
      report,
      file,
      format = "md"
    )
    
    lines <- readLines(
      file
    )
    
    testthat::expect_equal(
      lines[1],
      "# Test Report"
    )
    
  }
)


testthat::test_that(
  "write_report rejects invalid inputs",
  {
    
    report <- new_report(
      title = "Test Report",
      content = "Content"
    )
    
    testthat::expect_error(
      write_report(
        list(),
        tempfile(
          fileext = ".txt"
        )
      )
    )
    
    testthat::expect_error(
      write_report(
        report,
        NA_character_
      )
    )
    
    testthat::expect_error(
      write_report(
        report,
        tempfile(
          fileext = ".pdf"
        )
      )
    )
    
  }
)


testthat::test_that(
  "print.cpBarcodeReport prints a formatted report",
  {
    
    report <- new_report(
      title = "Test Report",
      content = c(
        "Line 1",
        "Line 2"
      )
    )
    
    output <- capture.output(
      print(
        report
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Test Report",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Line 1",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Line 2",
          output
        )
      )
    )
    
  }
)


testthat::test_that(
  "print.cpBarcodeReport rejects non-report objects",
  {
    
    not_report <- list(
      title = "Not a report",
      content = "Content"
    )
    
    class(not_report) <- "cpBarcodeReportInvalid"
    
    testthat::expect_error(
      print.cpBarcodeReport(
        not_report
      )
    )
    
  }
)