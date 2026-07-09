make_test_cpGenome <- function(){
  
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
    strand = "+",
    location_type = "simple",
    location = "1..8",
    gene = "TEST",
    stringsAsFactors = FALSE
  )
  
  new_cpGenome(
    metadata = metadata,
    sequence = sequence,
    feature_table = feature_table
  )
  
}


testthat::test_that(
  "new_cpGenome constructs a cpGenome object",
  {
    
    genome <- make_test_cpGenome()
    
    testthat::expect_true(
      inherits(
        genome,
        "cpGenome"
      )
    )
    
    testthat::expect_true(
      is.list(genome$metadata)
    )
    
    testthat::expect_equal(
      genome$sequence,
      "ATGCGTAA"
    )
    
    testthat::expect_true(
      is.data.frame(genome$feature_table)
    )
    
  }
)


testthat::test_that(
  "new_cpGenome rejects invalid constructor inputs",
  {
    
    valid_genome <- make_test_cpGenome()
    
    testthat::expect_error(
      new_cpGenome(
        metadata = "not metadata",
        sequence = valid_genome$sequence,
        feature_table = valid_genome$feature_table
      )
    )
    
    testthat::expect_error(
      new_cpGenome(
        metadata = valid_genome$metadata,
        sequence = c("ATGC", "GTAA"),
        feature_table = valid_genome$feature_table
      )
    )
    
    testthat::expect_error(
      new_cpGenome(
        metadata = valid_genome$metadata,
        sequence = valid_genome$sequence,
        feature_table = list()
      )
    )
    
  }
)


testthat::test_that(
  "validate_cpGenome rejects invalid cpGenome objects",
  {
    
    genome <- make_test_cpGenome()
    
    invalid_sequence <- genome
    invalid_sequence$sequence <- "ATGCXYZ"
    
    testthat::expect_error(
      validate_cpGenome(
        invalid_sequence
      )
    )
    
    empty_sequence <- genome
    empty_sequence$sequence <- ""
    
    testthat::expect_error(
      validate_cpGenome(
        empty_sequence
      )
    )
    
    invalid_features <- genome
    invalid_features$feature_table <- list()
    
    testthat::expect_error(
      validate_cpGenome(
        invalid_features
      )
    )
    
  }
)


testthat::test_that(
  "assert_cpGenome accepts cpGenome objects and rejects others",
  {
    
    genome <- make_test_cpGenome()
    
    testthat::expect_identical(
      assert_cpGenome(
        genome
      ),
      genome
    )
    
    testthat::expect_error(
      assert_cpGenome(
        list()
      )
    )
    
  }
)


testthat::test_that(
  "cpGenome accessors return expected values",
  {
    
    genome <- make_test_cpGenome()
    
    testthat::expect_equal(
      get_metadata(genome),
      genome$metadata
    )
    
    testthat::expect_equal(
      get_sequence(genome),
      "ATGCGTAA"
    )
    
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
  "print.cpGenome prints a concise genome summary",
  {
    
    genome <- make_test_cpGenome()
    
    output <- capture.output(
      print(
        genome
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "cpGenome Object",
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


genome_test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "read_metadata extracts genome metadata from a GenBank record",
  {
    
    gb_lines <- read_genbank(
      genome_test_file
    )
    
    metadata <- read_metadata(
      gb_lines
    )
    
    testthat::expect_true(
      is.list(metadata)
    )
    
    testthat::expect_true(
      all(
        c(
          "accession",
          "version",
          "organism",
          "definition",
          "genome_length"
        ) %in% names(metadata)
      )
    )
    
    testthat::expect_equal(
      metadata$accession,
      "NC_000932"
    )
    
    testthat::expect_equal(
      metadata$organism,
      "Arabidopsis thaliana"
    )
    
    testthat::expect_true(
      is.integer(metadata$genome_length)
    )
    
    testthat::expect_gt(
      metadata$genome_length,
      100000L
    )
    
  }
)


testthat::test_that(
  "read_cpGenome reads a GenBank file into a cpGenome object",
  {
    
    genome <- read_cpGenome(
      genome_test_file
    )
    
    testthat::expect_true(
      inherits(
        genome,
        "cpGenome"
      )
    )
    
    testthat::expect_equal(
      get_accession(genome),
      "NC_000932"
    )
    
    testthat::expect_equal(
      get_organism(genome),
      "Arabidopsis thaliana"
    )
    
    testthat::expect_true(
      is_single_string(
        get_sequence(genome)
      )
    )
    
    testthat::expect_gt(
      nchar(
        get_sequence(genome)
      ),
      100000L
    )
    
    testthat::expect_true(
      is.data.frame(
        get_feature_table(genome)
      )
    )
    
    testthat::expect_gt(
      nrow(
        get_feature_table(genome)
      ),
      0
    )
    
  }
)