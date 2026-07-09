pipeline_test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "extract_barcode extracts a single barcode from a GenBank file",
  {
    
    barcode <- extract_barcode(
      file = pipeline_test_file,
      marker = "RBCL"
    )
    
    testthat::expect_true(
      inherits(
        barcode,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_marker(barcode),
      MARKER_DICTIONARY$RBCL$marker
    )
    
    testthat::expect_equal(
      get_accession(barcode),
      "NC_000932"
    )
    
    testthat::expect_equal(
      get_organism(barcode),
      "Arabidopsis thaliana"
    )
    
    testthat::expect_gt(
      nchar(
        get_sequence(barcode)
      ),
      1000L
    )
    
  }
)


testthat::test_that(
  "extract_barcode writes output when requested",
  {
    
    output_file <- tempfile(
      fileext = ".fasta"
    )
    
    barcode <- extract_barcode(
      file = pipeline_test_file,
      marker = "RBCL",
      output = output_file,
      format = "fasta"
    )
    
    testthat::expect_true(
      inherits(
        barcode,
        "cpBarcode"
      )
    )
    
    testthat::expect_true(
      file.exists(output_file)
    )
    
    lines <- readLines(
      output_file
    )
    
    testthat::expect_true(
      length(lines) > 1
    )
    
    testthat::expect_true(
      grepl(
        "^>",
        lines[1]
      )
    )
    
  }
)


testthat::test_that(
  "extract_barcode rejects invalid marker input",
  {
    
    testthat::expect_error(
      extract_barcode(
        file = pipeline_test_file,
        marker = c("RBCL", "MATK")
      )
    )
    
    testthat::expect_error(
      extract_barcode(
        file = pipeline_test_file,
        marker = NA_character_
      )
    )
    
    testthat::expect_error(
      extract_barcode(
        file = pipeline_test_file,
        marker = "UNKNOWN_MARKER"
      )
    )
    
  }
)


testthat::test_that(
  "extract_barcodes extracts multiple barcode markers",
  {
    
    barcodes <- extract_barcodes(
      file = pipeline_test_file,
      markers = c(
        "RBCL",
        "MATK",
        "PSBA_TRNH"
      )
    )
    
    testthat::expect_true(
      is.list(barcodes)
    )
    
    testthat::expect_equal(
      names(barcodes),
      c(
        "RBCL",
        "MATK",
        "PSBA_TRNH"
      )
    )
    
    testthat::expect_equal(
      length(barcodes),
      3L
    )
    
    testthat::expect_true(
      all(
        vapply(
          barcodes,
          inherits,
          logical(1),
          what = "cpBarcode"
        )
      )
    )
    
    testthat::expect_equal(
      get_feature_type(
        barcodes$PSBA_TRNH
      ),
      FEATURE_INTERGENIC
    )
    
  }
)


testthat::test_that(
  "extract_barcodes writes output files when requested",
  {
    
    output_dir <- tempfile()
    
    barcodes <- extract_barcodes(
      file = pipeline_test_file,
      markers = c(
        "RBCL",
        "MATK"
      ),
      output_dir = output_dir,
      format = "fasta"
    )
    
    testthat::expect_true(
      is.list(barcodes)
    )
    
    testthat::expect_true(
      dir.exists(output_dir)
    )
    
    testthat::expect_true(
      file.exists(
        file.path(
          output_dir,
          "RBCL.fasta"
        )
      )
    )
    
    testthat::expect_true(
      file.exists(
        file.path(
          output_dir,
          "MATK.fasta"
        )
      )
    )
    
  }
)


testthat::test_that(
  "extract_barcodes rejects invalid inputs",
  {
    
    testthat::expect_error(
      extract_barcodes(
        file = pipeline_test_file,
        markers = character()
      )
    )
    
    testthat::expect_error(
      extract_barcodes(
        file = pipeline_test_file,
        markers = c(
          "RBCL",
          NA_character_
        )
      )
    )
    
    testthat::expect_error(
      extract_barcodes(
        file = pipeline_test_file,
        markers = "RBCL",
        output_dir = c(
          "a",
          "b"
        )
      )
    )
    
    testthat::expect_error(
      extract_barcodes(
        file = pipeline_test_file,
        markers = "UNKNOWN_MARKER"
      )
    )
    
  }
)