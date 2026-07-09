make_test_summary_barcode <- function(
    marker = MARKER_DICTIONARY[[1]]$marker,
    accession = "TEST001",
    organism = "Test organism",
    sequence = "AATTCCGG"
){
  
  new_cpBarcode(
    marker = marker,
    accession = accession,
    version = paste0(accession, ".1"),
    organism = organism,
    feature_type = FEATURE_GENE,
    sequence = sequence,
    coordinates = list(
      start = 1L,
      end = nchar(sequence)
    ),
    strand = STRAND_FORWARD,
    metadata = list()
  )
  
}


testthat::test_that(
  "summary.cpBarcode returns a cpBarcodeSummary object",
  {
    
    barcode <- make_test_summary_barcode()
    
    barcode_summary <- summary(
      barcode
    )
    
    testthat::expect_true(
      inherits(
        barcode_summary,
        "cpBarcodeSummary"
      )
    )
    
    testthat::expect_equal(
      barcode_summary$marker,
      get_marker(barcode)
    )
    
    testthat::expect_equal(
      barcode_summary$accession,
      "TEST001"
    )
    
    testthat::expect_equal(
      barcode_summary$organism,
      "Test organism"
    )
    
    testthat::expect_equal(
      barcode_summary$feature_type,
      FEATURE_GENE
    )
    
    testthat::expect_equal(
      barcode_summary$coordinates,
      list(
        start = 1L,
        end = 8L
      )
    )
    
    testthat::expect_equal(
      barcode_summary$strand,
      STRAND_FORWARD
    )
    
    testthat::expect_equal(
      barcode_summary$length,
      8L
    )
    
    testthat::expect_equal(
      barcode_summary$gc_content,
      50
    )
    
    testthat::expect_equal(
      barcode_summary$at_content,
      50
    )
    
    testthat::expect_equal(
      barcode_summary$ambiguous,
      0L
    )
    
    testthat::expect_equal(
      barcode_summary$ambiguous_percentage,
      0
    )
    
    testthat::expect_equal(
      barcode_summary$longest_homopolymer,
      2L
    )
    
  }
)


testthat::test_that(
  "print.cpBarcodeSummary prints a concise barcode summary",
  {
    
    barcode <- make_test_summary_barcode()
    
    barcode_summary <- summary(
      barcode
    )
    
    output <- capture.output(
      print(
        barcode_summary
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Barcode Summary",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Test organism",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "TEST001",
          output
        )
      )
    )
    
  }
)


testthat::test_that(
  "summarize_barcodes summarizes a collection of barcodes",
  {
    
    barcode1 <- make_test_summary_barcode(
      marker = MARKER_DICTIONARY[[1]]$marker,
      accession = "TEST001",
      organism = "Organism one",
      sequence = "AATTCCGG"
    )
    
    barcode2 <- make_test_summary_barcode(
      marker = MARKER_DICTIONARY[[2]]$marker,
      accession = "TEST002",
      organism = "Organism two",
      sequence = "AATTCCGGNN"
    )
    
    report <- summarize_barcodes(
      list(
        barcode1,
        barcode2
      )
    )
    
    testthat::expect_true(
      inherits(
        report,
        "cpBarcodeCollectionSummary"
      )
    )
    
    testthat::expect_equal(
      report$total_barcodes,
      2L
    )
    
    testthat::expect_equal(
      report$organism_count,
      2L
    )
    
    testthat::expect_equal(
      report$mean_length,
      9
    )
    
    testthat::expect_equal(
      report$min_length,
      8
    )
    
    testthat::expect_equal(
      report$max_length,
      10
    )
    
    testthat::expect_equal(
      report$total_ambiguous,
      2
    )
    
    testthat::expect_true(
      is.table(
        report$marker_frequency
      )
    )
    
  }
)


testthat::test_that(
  "summarize_barcodes rejects invalid inputs",
  {
    
    testthat::expect_error(
      summarize_barcodes(
        "not a list"
      )
    )
    
    testthat::expect_error(
      summarize_barcodes(
        list()
      )
    )
    
    testthat::expect_error(
      summarize_barcodes(
        list(
          list()
        )
      )
    )
    
  }
)


testthat::test_that(
  "print.cpBarcodeCollectionSummary prints a concise collection summary",
  {
    
    barcode1 <- make_test_summary_barcode(
      accession = "TEST001",
      organism = "Organism one",
      sequence = "AATTCCGG"
    )
    
    barcode2 <- make_test_summary_barcode(
      accession = "TEST002",
      organism = "Organism two",
      sequence = "AATTCCGGNN"
    )
    
    report <- summarize_barcodes(
      list(
        barcode1,
        barcode2
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
          "Barcode Collection Summary",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Total barcodes",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Marker frequencies",
          output
        )
      )
    )
    
  }
)