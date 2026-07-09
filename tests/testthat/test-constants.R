testthat::test_that(
  "package information constants are defined",
  {
    
    testthat::expect_equal(
      PACKAGE_NAME,
      "cpDNAbarcodes"
    )
    
    testthat::expect_equal(
      PACKAGE_VERSION,
      "1.0.0"
    )
    
  }
)


testthat::test_that(
  "GenBank section header constants are character patterns",
  {
    
    constants <- c(
      LOCUS_START,
      DEFINITION_START,
      ACCESSION_START,
      VERSION_START,
      ORGANISM_START,
      GENBANK_FEATURE_START,
      ORIGIN_START,
      GENBANK_END
    )
    
    testthat::expect_true(
      all(
        vapply(
          constants,
          is_single_string,
          logical(1)
        )
      )
    )
    
    testthat::expect_true(
      grepl(
        LOCUS_START,
        "LOCUS       NC_000932"
      )
    )
    
    testthat::expect_true(
      grepl(
        ACCESSION_START,
        "ACCESSION   NC_000932"
      )
    )
    
    testthat::expect_true(
      grepl(
        ORGANISM_START,
        "  ORGANISM  Arabidopsis thaliana"
      )
    )
    
    testthat::expect_true(
      grepl(
        ORIGIN_START,
        "ORIGIN"
      )
    )
    
    testthat::expect_true(
      grepl(
        GENBANK_END,
        "//"
      )
    )
    
  }
)


testthat::test_that(
  "feature type constants are defined",
  {
    
    testthat::expect_equal(
      FEATURE_GENE,
      "gene"
    )
    
    testthat::expect_equal(
      FEATURE_CDS,
      "CDS"
    )
    
    testthat::expect_equal(
      FEATURE_TRNA,
      "tRNA"
    )
    
    testthat::expect_equal(
      FEATURE_RRNA,
      "rRNA"
    )
    
    testthat::expect_equal(
      FEATURE_INTERGENIC,
      "intergenic"
    )
    
  }
)


testthat::test_that(
  "marker class constants are defined",
  {
    
    testthat::expect_equal(
      MARKER_CLASS_GENE,
      "GENE"
    )
    
    testthat::expect_equal(
      MARKER_CLASS_INTERGENIC,
      "INTERGENIC"
    )
    
    testthat::expect_equal(
      MARKER_CLASS_NUCLEAR,
      "NUCLEAR"
    )
    
  }
)


testthat::test_that(
  "strand constants are defined",
  {
    
    testthat::expect_equal(
      STRAND_FORWARD,
      "+"
    )
    
    testthat::expect_equal(
      STRAND_REVERSE,
      "-"
    )
    
    testthat::expect_true(
      all(
        c(
          STRAND_FORWARD,
          STRAND_REVERSE
        ) %in% c(
          "+",
          "-"
        )
      )
    )
    
  }
)


testthat::test_that(
  "location type constants are defined",
  {
    
    testthat::expect_equal(
      LOCATION_SIMPLE,
      "simple"
    )
    
    testthat::expect_equal(
      LOCATION_JOIN,
      "join"
    )
    
    testthat::expect_equal(
      LOCATION_ORDER,
      "order"
    )
    
  }
)


testthat::test_that(
  "DNA alphabet constants are defined",
  {
    
    testthat::expect_equal(
      DNA_BASES,
      c(
        "A",
        "C",
        "G",
        "T"
      )
    )
    
    testthat::expect_true(
      all(
        DNA_BASES %in% IUPAC_BASES
      )
    )
    
    testthat::expect_true(
      "N" %in% IUPAC_BASES
    )
    
    testthat::expect_equal(
      length(
        unique(
          IUPAC_BASES
        )
      ),
      length(
        IUPAC_BASES
      )
    )
    
  }
)


testthat::test_that(
  "export default constants are defined",
  {
    
    testthat::expect_equal(
      DEFAULT_FASTA_WIDTH,
      80
    )
    
    testthat::expect_equal(
      DEFAULT_FASTA_EXTENSION,
      ".fasta"
    )
    
    testthat::expect_true(
      DEFAULT_FASTA_WIDTH > 0
    )
    
  }
)