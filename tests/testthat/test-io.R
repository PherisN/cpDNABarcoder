test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "read_genbank reads a valid GenBank file",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    testthat::expect_true(
      is.character(lines)
    )
    
    testthat::expect_gt(
      length(lines),
      0
    )
    
  }
)


testthat::test_that(
  "extract_accession returns the accession number",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    accession <- extract_accession(
      lines
    )
    
    testthat::expect_true(
      is_single_string(accession)
    )
    
    testthat::expect_equal(
      accession,
      "NC_000932"
    )
    
  }
)


testthat::test_that(
  "extract_organism returns the organism name",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    organism <- extract_organism(
      lines
    )
    
    testthat::expect_true(
      is_single_string(organism)
    )
    
    testthat::expect_equal(
      organism,
      "Arabidopsis thaliana"
    )
    
  }
)


testthat::test_that(
  "extract_features returns the FEATURES section",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    features <- extract_features(
      lines
    )
    
    testthat::expect_true(
      is.character(features)
    )
    
    testthat::expect_gt(
      length(features),
      0
    )
    
  }
)


testthat::test_that(
  "read_sequence returns a valid DNA sequence",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    sequence <- read_sequence(
      lines
    )
    
    testthat::expect_true(
      is_single_string(sequence)
    )
    
    testthat::expect_true(
      is_valid_dna(sequence)
    )
    
    testthat::expect_true(
      nchar(sequence) > 100000
    )
    
  }
)


testthat::test_that(
  "read_genbank rejects invalid file arguments",
  {
    
    testthat::expect_error(
      read_genbank(1)
    )
    
    testthat::expect_error(
      read_genbank("this_file_does_not_exist.gb")
    )
    
    testthat::expect_error(
      read_genbank(project_root)
    )
    
  }
)


testthat::test_that(
  "extract_accession rejects invalid input",
  {
    
    testthat::expect_error(
      extract_accession(1)
    )
    
    testthat::expect_error(
      extract_accession(character())
    )
    
  }
)


testthat::test_that(
  "extract_organism rejects invalid input",
  {
    
    testthat::expect_error(
      extract_organism(1)
    )
    
    testthat::expect_error(
      extract_organism(character())
    )
    
  }
)


testthat::test_that(
  "extract_features rejects invalid input",
  {
    
    testthat::expect_error(
      extract_features(1)
    )
    
    testthat::expect_error(
      extract_features(character())
    )
    
  }
)


testthat::test_that(
  "read_sequence rejects invalid input",
  {
    
    testthat::expect_error(
      read_sequence(1)
    )
    
    testthat::expect_error(
      read_sequence(character())
    )
    
  }
)