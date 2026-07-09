make_test_qc_barcode <- function(){
  
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
  "is_valid_dna recognises valid and invalid DNA sequences",
  {
    
    testthat::expect_true(
      is_valid_dna(
        "ATGCCGTAN"
      )
    )
    
    testthat::expect_true(
      is_valid_dna(
        "atgccgtan"
      )
    )
    
    testthat::expect_false(
      is_valid_dna(
        "ATGCX"
      )
    )
    
    testthat::expect_false(
      is_valid_dna(
        c("ATGC", "ATGC")
      )
    )
    
    testthat::expect_false(
      is_valid_dna(
        ""
      )
    )
    
  }
)


testthat::test_that(
  "count_ambiguous_bases counts non-ACGT symbols",
  {
    
    testthat::expect_equal(
      count_ambiguous_bases(
        "AATTCCGGNN"
      ),
      2L
    )
    
    testthat::expect_equal(
      count_ambiguous_bases(
        "AATTCCGG"
      ),
      0L
    )
    
    testthat::expect_error(
      count_ambiguous_bases(
        "ATGCX"
      )
    )
    
  }
)


testthat::test_that(
  "nucleotide_composition counts nucleotide symbols",
  {
    
    composition <- nucleotide_composition(
      "AATTCCGGNN"
    )
    
    testthat::expect_equal(
      as.integer(
        composition["A"]
      ),
      2L
    )
    
    testthat::expect_equal(
      as.integer(
        composition["T"]
      ),
      2L
    )
    
    testthat::expect_equal(
      as.integer(
        composition["C"]
      ),
      2L
    )
    
    testthat::expect_equal(
      as.integer(
        composition["G"]
      ),
      2L
    )
    
    testthat::expect_equal(
      as.integer(
        composition["N"]
      ),
      2L
    )
    
  }
)


testthat::test_that(
  "gc_content and at_content calculate nucleotide percentages",
  {
    
    sequence <- "AATTCCGGNN"
    
    testthat::expect_equal(
      gc_content(sequence),
      40
    )
    
    testthat::expect_equal(
      at_content(sequence),
      40
    )
    
    testthat::expect_error(
      gc_content(
        "ATGCX"
      )
    )
    
    testthat::expect_error(
      at_content(
        "ATGCX"
      )
    )
    
  }
)


testthat::test_that(
  "ambiguous_percentage calculates ambiguous base percentage",
  {
    
    testthat::expect_equal(
      ambiguous_percentage(
        "AATTCCGGNN"
      ),
      20
    )
    
    testthat::expect_equal(
      ambiguous_percentage(
        "AATTCCGG"
      ),
      0
    )
    
    testthat::expect_error(
      ambiguous_percentage(
        "ATGCX"
      )
    )
    
  }
)


testthat::test_that(
  "longest_homopolymer returns the longest repeated base run",
  {
    
    testthat::expect_equal(
      longest_homopolymer(
        "AAATTTCCG"
      ),
      3L
    )
    
    testthat::expect_equal(
      longest_homopolymer(
        "ATGC"
      ),
      1L
    )
    
    testthat::expect_error(
      longest_homopolymer(
        "ATGCX"
      )
    )
    
  }
)


testthat::test_that(
  "sequence_length returns DNA sequence length",
  {
    
    testthat::expect_equal(
      sequence_length(
        "ATGCNN"
      ),
      6L
    )
    
    testthat::expect_error(
      sequence_length(
        "ATGCX"
      )
    )
    
  }
)


testthat::test_that(
  "qc_barcode returns a cpBarcodeQC report",
  {
    
    barcode <- make_test_qc_barcode()
    
    report <- qc_barcode(
      barcode
    )
    
    testthat::expect_true(
      inherits(
        report,
        "cpBarcodeQC"
      )
    )
    
    testthat::expect_equal(
      report$marker,
      get_marker(barcode)
    )
    
    testthat::expect_equal(
      report$accession,
      "TEST001"
    )
    
    testthat::expect_true(
      report$valid
    )
    
    testthat::expect_equal(
      report$length,
      10L
    )
    
    testthat::expect_equal(
      report$gc_content,
      40
    )
    
    testthat::expect_equal(
      report$at_content,
      40
    )
    
    testthat::expect_equal(
      report$ambiguous,
      2L
    )
    
    testthat::expect_equal(
      report$ambiguous_percentage,
      20
    )
    
    testthat::expect_equal(
      report$longest_homopolymer,
      2L
    )
    
  }
)


testthat::test_that(
  "qc_barcode records warnings for poor-quality sequences",
  {
    
    barcode <- make_test_qc_barcode()
    
    report <- qc_barcode(
      barcode
    )
    
    testthat::expect_true(
      "High proportion of ambiguous bases." %in%
        report$warnings
    )
    
  }
)


testthat::test_that(
  "assert_cpBarcodeQC accepts QC reports and rejects others",
  {
    
    barcode <- make_test_qc_barcode()
    
    report <- qc_barcode(
      barcode
    )
    
    testthat::expect_identical(
      assert_cpBarcodeQC(
        report
      ),
      report
    )
    
    testthat::expect_error(
      assert_cpBarcodeQC(
        list()
      )
    )
    
  }
)


testthat::test_that(
  "print.cpBarcodeQC prints a concise QC summary",
  {
    
    barcode <- make_test_qc_barcode()
    
    report <- qc_barcode(
      barcode
    )
    
    output <- capture.output(
      print(
        report
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Barcode QC Report",
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
    
    testthat::expect_true(
      any(
        grepl(
          "High proportion of ambiguous bases",
          output
        )
      )
    )
    
  }
)