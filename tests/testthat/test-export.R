make_test_export_barcode <- function(){
  
  valid_marker <- MARKER_DICTIONARY[[1]]$marker
  
  new_cpBarcode(
    marker = valid_marker,
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    feature_type = FEATURE_GENE,
    sequence = "AATTCCGGNN",
    coordinates = list(
      start = 10L,
      end = 19L
    ),
    strand = STRAND_FORWARD,
    metadata = list()
  )
  
}


testthat::test_that(
  "wrap_fasta_sequence wraps DNA sequences to the requested width",
  {
    
    wrapped <- wrap_fasta_sequence(
      "AATTCCGGNN",
      width = 4
    )
    
    testthat::expect_equal(
      wrapped,
      c(
        "AATT",
        "CCGG",
        "NN"
      )
    )
    
  }
)


testthat::test_that(
  "wrap_fasta_sequence rejects invalid inputs",
  {
    
    testthat::expect_error(
      wrap_fasta_sequence(
        "ATGCX",
        width = 4
      )
    )
    
    testthat::expect_error(
      wrap_fasta_sequence(
        "",
        width = 4
      )
    )
    
    testthat::expect_error(
      wrap_fasta_sequence(
        "ATGC",
        width = 0
      )
    )
    
  }
)


testthat::test_that(
  "build_fasta_header builds a standardized FASTA header",
  {
    
    barcode <- make_test_export_barcode()
    
    header <- build_fasta_header(
      barcode
    )
    
    testthat::expect_equal(
      header,
      paste0(
        ">",
        get_accession(barcode),
        "|",
        get_marker(barcode),
        "|Test organism|10-19"
      )
    )
    
  }
)


testthat::test_that(
  "as_fasta converts a barcode to FASTA text",
  {
    
    barcode <- make_test_export_barcode()
    
    fasta <- as_fasta(
      barcode,
      width = 5
    )
    
    testthat::expect_true(
      is.character(fasta)
    )
    
    testthat::expect_equal(
      fasta[1],
      build_fasta_header(barcode)
    )
    
    testthat::expect_equal(
      fasta[-1],
      c(
        "AATTC",
        "CGGNN"
      )
    )
    
  }
)


testthat::test_that(
  "write_fasta writes a FASTA file",
  {
    
    barcode <- make_test_export_barcode()
    
    file <- tempfile(
      fileext = ".fasta"
    )
    
    result <- write_fasta(
      barcode,
      file,
      width = 5
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
      lines,
      as_fasta(
        barcode,
        width = 5
      )
    )
    
  }
)


testthat::test_that(
  "as_dataframe converts a barcode to a one-row data frame",
  {
    
    barcode <- make_test_export_barcode()
    
    table <- as_dataframe(
      barcode
    )
    
    testthat::expect_true(
      is.data.frame(table)
    )
    
    testthat::expect_equal(
      nrow(table),
      1L
    )
    
    testthat::expect_equal(
      table$marker,
      get_marker(barcode)
    )
    
    testthat::expect_equal(
      table$accession,
      "TEST001"
    )
    
    testthat::expect_equal(
      table$version,
      "TEST001.1"
    )
    
    testthat::expect_equal(
      table$organism,
      "Test organism"
    )
    
    testthat::expect_equal(
      table$feature_type,
      FEATURE_GENE
    )
    
    testthat::expect_equal(
      table$start,
      10L
    )
    
    testthat::expect_equal(
      table$end,
      19L
    )
    
    testthat::expect_equal(
      table$strand,
      STRAND_FORWARD
    )
    
    testthat::expect_equal(
      table$length,
      10L
    )
    
    testthat::expect_equal(
      table$sequence,
      "AATTCCGGNN"
    )
    
  }
)


testthat::test_that(
  "write_csv writes a CSV file",
  {
    
    barcode <- make_test_export_barcode()
    
    file <- tempfile(
      fileext = ".csv"
    )
    
    result <- write_csv(
      barcode,
      file
    )
    
    testthat::expect_equal(
      result,
      file
    )
    
    testthat::expect_true(
      file.exists(file)
    )
    
    table <- read.csv(
      file,
      stringsAsFactors = FALSE
    )
    
    testthat::expect_equal(
      nrow(table),
      1L
    )
    
    testthat::expect_equal(
      table$accession,
      "TEST001"
    )
    
    testthat::expect_equal(
      table$sequence,
      "AATTCCGGNN"
    )
    
  }
)


testthat::test_that(
  "write_barcode dispatches supported export formats",
  {
    
    barcode <- make_test_export_barcode()
    
    fasta_file <- tempfile(
      fileext = ".fasta"
    )
    
    csv_file <- tempfile(
      fileext = ".csv"
    )
    
    write_barcode(
      barcode,
      fasta_file,
      format = "fasta",
      width = 5
    )
    
    write_barcode(
      barcode,
      csv_file,
      format = "csv"
    )
    
    testthat::expect_true(
      file.exists(fasta_file)
    )
    
    testthat::expect_true(
      file.exists(csv_file)
    )
    
  }
)


testthat::test_that(
  "write_barcode normalizes format names",
  {
    
    barcode <- make_test_export_barcode()
    
    fasta_file <- tempfile(
      fileext = ".fa"
    )
    
    write_barcode(
      barcode,
      fasta_file,
      format = ".FASTA",
      width = 5
    )
    
    testthat::expect_true(
      file.exists(fasta_file)
    )
    
  }
)


testthat::test_that(
  "write_barcode rejects unsupported formats",
  {
    
    barcode <- make_test_export_barcode()
    
    file <- tempfile(
      fileext = ".txt"
    )
    
    testthat::expect_error(
      write_barcode(
        barcode,
        file,
        format = "txt"
      )
    )
    
  }
)