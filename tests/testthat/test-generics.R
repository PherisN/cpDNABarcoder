make_test_generics_genome <- function(){
  
  metadata <- list(
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    definition = "Test chloroplast genome",
    genome_length = 8L
  )
  
  sequence <- "ATGCGTAA"
  
  feature_table <- data.frame(
    type = "gene",
    start = 1L,
    end = 8L,
    strand = STRAND_FORWARD,
    location_type = LOCATION_SIMPLE,
    location = "1..8",
    gene = "rbcL",
    stringsAsFactors = FALSE
  )
  
  new_cpGenome(
    metadata = metadata,
    sequence = sequence,
    feature_table = feature_table
  )
  
}


make_test_generics_barcode <- function(){
  
  valid_marker <- MARKER_DICTIONARY[[1]]$marker
  
  new_cpBarcode(
    marker = valid_marker,
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    feature_type = FEATURE_GENE,
    sequence = "ATGCGTAA",
    coordinates = list(
      start = 1L,
      end = 8L
    ),
    strand = STRAND_FORWARD,
    metadata = list(
      source = "unit test"
    )
  )
  
}


testthat::test_that(
  "shared generics dispatch for cpGenome objects",
  {
    
    genome <- make_test_generics_genome()
    
    testthat::expect_equal(
      get_metadata(genome),
      genome$metadata
    )
    
    testthat::expect_equal(
      get_sequence(genome),
      "ATGCGTAA"
    )
    
  }
)


testthat::test_that(
  "genome generics dispatch for cpGenome objects",
  {
    
    genome <- make_test_generics_genome()
    
    testthat::expect_equal(
      get_feature_table(genome),
      genome$feature_table
    )
    
    testthat::expect_equal(
      get_accession(genome),
      "TEST001"
    )
    
    testthat::expect_equal(
      get_version(genome),
      "TEST001.1"
    )
    
    testthat::expect_equal(
      get_definition(genome),
      "Test chloroplast genome"
    )
    
    testthat::expect_equal(
      get_organism(genome),
      "Test organism"
    )
    
    testthat::expect_equal(
      get_genome_length(genome),
      8L
    )
    
  }
)


testthat::test_that(
  "shared generics dispatch for cpBarcode objects",
  {
    
    barcode <- make_test_generics_barcode()
    
    testthat::expect_equal(
      get_metadata(barcode),
      list(
        source = "unit test"
      )
    )
    
    testthat::expect_equal(
      get_sequence(barcode),
      "ATGCGTAA"
    )
    
  }
)


testthat::test_that(
  "barcode generics dispatch for cpBarcode objects",
  {
    
    barcode <- make_test_generics_barcode()
    
    testthat::expect_equal(
      get_marker(barcode),
      barcode$marker
    )
    
    testthat::expect_equal(
      get_coordinates(barcode),
      list(
        start = 1L,
        end = 8L
      )
    )
    
    testthat::expect_equal(
      get_strand(barcode),
      STRAND_FORWARD
    )
    
    testthat::expect_equal(
      get_feature_type(barcode),
      FEATURE_GENE
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
    
  }
)


testthat::test_that(
  "generics reject unsupported object classes",
  {
    
    unsupported <- list()
    
    testthat::expect_error(
      get_metadata(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_sequence(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_accession(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_marker(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_coordinates(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_strand(
        unsupported
      )
    )
    
    testthat::expect_error(
      get_feature_type(
        unsupported
      )
    )
    
  }
)