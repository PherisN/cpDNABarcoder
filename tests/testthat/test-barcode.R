make_test_cpBarcode <- function(){
  
  valid_marker <- MARKER_DICTIONARY[[1]]$marker
  
  new_cpBarcode(
    marker = valid_marker,
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    feature_type = FEATURE_GENE,
    sequence = "ATGCGTAA",
    coordinates = list(
      start = 10L,
      end = 17L
    ),
    strand = STRAND_FORWARD,
    metadata = list(
      source = "unit test"
    )
  )
  
}


testthat::test_that(
  "new_cpBarcode constructs a cpBarcode object",
  {
    
    barcode <- make_test_cpBarcode()
    
    testthat::expect_true(
      inherits(
        barcode,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      barcode$sequence,
      "ATGCGTAA"
    )
    
    testthat::expect_equal(
      barcode$length,
      8L
    )
    
    testthat::expect_equal(
      barcode$coordinates$start,
      10L
    )
    
    testthat::expect_equal(
      barcode$coordinates$end,
      17L
    )
    
    testthat::expect_equal(
      barcode$strand,
      STRAND_FORWARD
    )
    
  }
)


testthat::test_that(
  "new_cpBarcode rejects invalid constructor inputs",
  {
    
    valid_marker <- MARKER_DICTIONARY[[1]]$marker
    
    testthat::expect_error(
      new_cpBarcode(
        marker = valid_marker,
        accession = "TEST001",
        version = "TEST001.1",
        organism = "Test organism",
        feature_type = FEATURE_GENE,
        sequence = c("ATGC", "GTAA"),
        coordinates = list(
          start = 10L,
          end = 17L
        ),
        strand = STRAND_FORWARD
      )
    )
    
    testthat::expect_error(
      new_cpBarcode(
        marker = valid_marker,
        accession = "TEST001",
        version = "TEST001.1",
        organism = "Test organism",
        feature_type = FEATURE_GENE,
        sequence = "",
        coordinates = list(
          start = 10L,
          end = 17L
        ),
        strand = STRAND_FORWARD
      )
    )
    
    testthat::expect_error(
      new_cpBarcode(
        marker = valid_marker,
        accession = "TEST001",
        version = "TEST001.1",
        organism = "Test organism",
        feature_type = FEATURE_GENE,
        sequence = "ATGCGTAA",
        coordinates = list(
          start = NA_integer_,
          end = 17L
        ),
        strand = STRAND_FORWARD
      )
    )
    
    testthat::expect_error(
      new_cpBarcode(
        marker = valid_marker,
        accession = "TEST001",
        version = "TEST001.1",
        organism = "Test organism",
        feature_type = FEATURE_GENE,
        sequence = "ATGCGTAA",
        coordinates = list(
          start = 17L,
          end = 10L
        ),
        strand = STRAND_FORWARD
      )
    )
    
  }
)


testthat::test_that(
  "assert_cpBarcode accepts cpBarcode objects and rejects others",
  {
    
    barcode <- make_test_cpBarcode()
    
    testthat::expect_identical(
      assert_cpBarcode(
        barcode
      ),
      barcode
    )
    
    testthat::expect_error(
      assert_cpBarcode(
        list()
      )
    )
    
  }
)


testthat::test_that(
  "validate_barcode rejects invalid cpBarcode objects",
  {
    
    barcode <- make_test_cpBarcode()
    
    invalid_sequence <- barcode
    invalid_sequence$sequence <- "ATGCXYZ"
    
    testthat::expect_error(
      validate_barcode(
        invalid_sequence
      )
    )
    
    invalid_length <- barcode
    invalid_length$length <- 999L
    
    testthat::expect_error(
      validate_barcode(
        invalid_length
      )
    )
    
    invalid_strand <- barcode
    invalid_strand$strand <- "?"
    
    testthat::expect_error(
      validate_barcode(
        invalid_strand
      )
    )
    
    invalid_metadata <- barcode
    invalid_metadata$metadata <- "not metadata"
    
    testthat::expect_error(
      validate_barcode(
        invalid_metadata
      )
    )
    
    missing_field <- barcode
    missing_field$marker <- NULL
    
    testthat::expect_error(
      validate_barcode(
        missing_field
      )
    )
    
  }
)


testthat::test_that(
  "validate_barcode warns for unknown feature types and markers",
  {
    
    barcode <- make_test_cpBarcode()
    
    unknown_feature_type <- barcode
    unknown_feature_type$feature_type <- "unknown_feature"
    
    testthat::expect_warning(
      validate_barcode(
        unknown_feature_type
      )
    )
    
    unknown_marker <- barcode
    unknown_marker$marker <- "UNKNOWN_MARKER"
    
    testthat::expect_warning(
      validate_barcode(
        unknown_marker
      )
    )
    
  }
)


testthat::test_that(
  "cpBarcode accessors return expected values",
  {
    
    barcode <- make_test_cpBarcode()
    
    testthat::expect_equal(
      get_sequence(barcode),
      "ATGCGTAA"
    )
    
    testthat::expect_equal(
      get_coordinates(barcode),
      list(
        start = 10L,
        end = 17L
      )
    )
    
    testthat::expect_equal(
      get_metadata(barcode),
      list(
        source = "unit test"
      )
    )
    
    testthat::expect_equal(
      get_marker(barcode),
      barcode$marker
    )
    
    testthat::expect_equal(
      get_accession(barcode),
      "TEST001"
    )
    
    testthat::expect_equal(
      get_version(barcode),
      "TEST001.1"
    )
    
    testthat::expect_equal(
      get_organism(barcode),
      "Test organism"
    )
    
    testthat::expect_equal(
      get_feature_type(barcode),
      FEATURE_GENE
    )
    
    testthat::expect_equal(
      get_strand(barcode),
      STRAND_FORWARD
    )
    
  }
)


testthat::test_that(
  "print.cpBarcode prints a concise barcode summary",
  {
    
    barcode <- make_test_cpBarcode()
    
    output <- capture.output(
      print(
        barcode
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "cpBarcode Object",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          barcode$marker,
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