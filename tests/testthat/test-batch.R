batch_test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "batch_extract_barcode extracts one marker from multiple files",
  {
    
    results <- suppressMessages(
      batch_extract_barcode(
        files = c(
          batch_test_file,
          batch_test_file
        ),
        marker = "RBCL"
      )
    )
    
    testthat::expect_true(
      is.list(results)
    )
    
    testthat::expect_equal(
      length(results),
      2L
    )
    
    testthat::expect_true(
      all(
        vapply(
          results,
          inherits,
          logical(1),
          what = "cpBarcode"
        )
      )
    )
    
    testthat::expect_true(
      all(
        vapply(
          results,
          function(x) get_marker(x) == MARKER_DICTIONARY$RBCL$marker,
          logical(1)
        )
      )
    )
    
  }
)


testthat::test_that(
  "batch_extract_barcode returns NULL for failed files without stopping",
  {
    
    missing_file <- file.path(
      tempdir(),
      "missing_file.gb"
    )
    
    results <- NULL
    
    testthat::expect_warning(
      results <- suppressMessages(
        batch_extract_barcode(
          files = c(
            batch_test_file,
            missing_file
          ),
          marker = "RBCL"
        )
      )
    )
    
    testthat::expect_equal(
      length(results),
      2L
    )
    
    testthat::expect_true(
      inherits(
        results[[1]],
        "cpBarcode"
      )
    )
    
    testthat::expect_null(
      results[[2]]
    )
    
  }
)


testthat::test_that(
  "batch_extract_barcode writes output files when requested",
  {
    
    output_dir <- tempfile()
    
    results <- suppressMessages(
      batch_extract_barcode(
        files = batch_test_file,
        marker = "RBCL",
        output_dir = output_dir,
        format = "fasta"
      )
    )
    
    expected_file <- file.path(
      output_dir,
      paste0(
        tools::file_path_sans_ext(
          basename(batch_test_file)
        ),
        ".fasta"
      )
    )
    
    testthat::expect_true(
      inherits(
        results[[1]],
        "cpBarcode"
      )
    )
    
    testthat::expect_true(
      dir.exists(output_dir)
    )
    
    testthat::expect_true(
      file.exists(expected_file)
    )
    
    lines <- readLines(
      expected_file
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
  "batch_extract_barcode rejects invalid inputs",
  {
    
    testthat::expect_error(
      batch_extract_barcode(
        files = character(),
        marker = "RBCL"
      )
    )
    
    testthat::expect_error(
      batch_extract_barcode(
        files = c(
          batch_test_file,
          NA_character_
        ),
        marker = "RBCL"
      )
    )
    
    testthat::expect_error(
      batch_extract_barcode(
        files = batch_test_file,
        marker = c(
          "RBCL",
          "MATK"
        )
      )
    )
    
    testthat::expect_error(
      batch_extract_barcode(
        files = batch_test_file,
        marker = "RBCL",
        output_dir = c(
          "a",
          "b"
        )
      )
    )
    
    testthat::expect_error(
      batch_extract_barcode(
        files = batch_test_file,
        marker = "RBCL",
        format = c(
          "fasta",
          "csv"
        )
      )
    )
    
  }
)


testthat::test_that(
  "batch_extract_barcodes extracts multiple markers from files",
  {
    
    results <- suppressMessages(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = c(
          "RBCL",
          "MATK",
          "PSBA_TRNH"
        )
      )
    )
    
    testthat::expect_true(
      is.list(results)
    )
    
    testthat::expect_equal(
      length(results),
      1L
    )
    
    testthat::expect_true(
      is.list(
        results[[1]]
      )
    )
    
    testthat::expect_equal(
      names(
        results[[1]]
      ),
      c(
        "RBCL",
        "MATK",
        "PSBA_TRNH"
      )
    )
    
    testthat::expect_true(
      all(
        vapply(
          results[[1]],
          inherits,
          logical(1),
          what = "cpBarcode"
        )
      )
    )
    
  }
)


testthat::test_that(
  "batch_extract_barcodes returns NULL for failed files without stopping",
  {
    
    missing_file <- file.path(
      tempdir(),
      "missing_file.gb"
    )
    
    results <- NULL
    
    testthat::expect_warning(
      results <- suppressMessages(
        batch_extract_barcodes(
          files = c(
            batch_test_file,
            missing_file
          ),
          markers = c(
            "RBCL",
            "MATK"
          )
        )
      )
    )
    
    testthat::expect_equal(
      length(results),
      2L
    )
    
    testthat::expect_true(
      is.list(
        results[[1]]
      )
    )
    
    testthat::expect_null(
      results[[2]]
    )
    
  }
)


testthat::test_that(
  "batch_extract_barcodes writes output files when requested",
  {
    
    output_dir <- tempfile()
    
    results <- suppressMessages(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = c(
          "RBCL",
          "MATK"
        ),
        output_dir = output_dir,
        format = "fasta"
      )
    )
    
    testthat::expect_true(
      is.list(
        results[[1]]
      )
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
  "batch_extract_barcodes rejects invalid inputs",
  {
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = character(),
        markers = "RBCL"
      )
    )
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = c(
          batch_test_file,
          NA_character_
        ),
        markers = "RBCL"
      )
    )
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = character()
      )
    )
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = c(
          "RBCL",
          NA_character_
        )
      )
    )
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = "RBCL",
        output_dir = c(
          "a",
          "b"
        )
      )
    )
    
    testthat::expect_error(
      batch_extract_barcodes(
        files = batch_test_file,
        markers = "RBCL",
        format = c(
          "fasta",
          "csv"
        )
      )
    )
    
  }
)


test_that(
  "batch_extract_barcodes writes all genomes to combined marker FASTA files",
  {
    
    input_dir <- tempfile(
      "combined_batch_inputs_"
    )
    
    output_dir <- tempfile(
      "combined_batch_fasta_"
    )
    
    dir.create(
      input_dir
    )
    
    dir.create(
      output_dir
    )
    
    on.exit(
      unlink(
        c(
          input_dir,
          output_dir
        ),
        recursive = TRUE,
        force = TRUE
      ),
      add = TRUE
    )
    
    source_file <- test_genbank_file(
      "NC_000932.gb"
    )
    
    files <- file.path(
      input_dir,
      c(
        "genome_1.gb",
        "genome_2.gb"
      )
    )
    
    file.copy(
      from = source_file,
      to = files
    )
    
    results <- batch_extract_barcodes(
      files = files,
      markers = c(
        "RBCL",
        "MATK",
        "PSBA_TRNH"
      ),
      output_dir = output_dir,
      format = "fasta"
    )
    
    expect_length(
      results,
      2
    )
    
    output_files <- file.path(
      output_dir,
      c(
        "RBCL.fasta",
        "MATK.fasta",
        "PSBA_TRNH.fasta"
      )
    )
    
    expect_true(
      all(
        file.exists(
          output_files
        )
      )
    )
    
    sequence_counts <- vapply(
      output_files,
      function(file){
        
        sum(
          grepl(
            "^>",
            readLines(
              file,
              warn = FALSE
            )
          )
        )
        
      },
      integer(1)
    )
    
    expect_equal(
      unname(
        sequence_counts
      ),
      c(
        2L,
        2L,
        2L
      )
    )
    
  }
)
